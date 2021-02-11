//BinaryDump::BinaryDump(int _rank, int _nproc) : Dump_Strategy(_rank, _nproc)
//{
    //// empty
//}
#include "dump_strategy.h"

void BinaryDump::dump_fields(
        const char *fbase,
        int step,
        grid_t* grid,
        field_array_t* field_array,
        int ftag
        )
{
    char fname[256];
    FileIO fileIO;
    int dim[3];

    if( !fbase ) ERROR(( "Invalid filename" ));

    if( rank==0 ) MESSAGE(( "Dumping fields to \"%s\"", fbase ));

    if( ftag ) sprintf( fname, "%s.%li.%i", fbase, (long)step, rank );
    else       sprintf( fname, "%s.%i", fbase, rank );

    FileIOStatus status = fileIO.open(fname, io_write);
    if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));

    /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
    size_t nxout = grid->nx;
    size_t nyout = grid->ny;
    size_t nzout = grid->nz;
    float dxout = grid->dx;
    float dyout = grid->dy;
    float dzout = grid->dz;

    WRITE_HEADER_V0( dump_type::field_dump, -1, 0, fileIO, step , rank, nproc);

    dim[0] = grid->nx+2;
    dim[1] = grid->ny+2;
    dim[2] = grid->nz+2;
    WRITE_ARRAY_HEADER( field_array->f, 3, dim, fileIO );
    fileIO.write( field_array->f, dim[0]*dim[1]*dim[2] );
    if( fileIO.close() ) ERROR(( "File close failed on dump fields!!!" ));
}

void BinaryDump::dump_particles(
        const char *fbase,
        species_t* sp,
        grid_t* grid,
        int step,
        interpolator_array_t* interpolator_array,
        int ftag
        )
{
    char fname[256];
    FileIO fileIO;
    int dim[1], buf_start;
    static particle_t * ALIGNED(128) p_buf = NULL;

    // TODO: reconcile this with MAX_IO_CHUNK, and update Cmake option
    // description to explain what backends use it
# define PBUF_SIZE 32768 // 1MB of particles

    if( !sp ) ERROR(( "Invalid species name \"%s\".", sp->name ));

    if( !fbase ) ERROR(( "Invalid filename" ));

    if( !p_buf ) MALLOC_ALIGNED( p_buf, PBUF_SIZE, 128 );

    if( rank==0 )
        MESSAGE(("Dumping \"%s\" particles to \"%s\"",sp->name,fbase));

    if( ftag ) sprintf( fname, "%s.%li.%i", fbase, (long)step, rank );
    else       sprintf( fname, "%s.%i", fbase, rank );
    FileIOStatus status = fileIO.open(fname, io_write);
    if( status==fail ) ERROR(( "Could not open \"%s\"", fname ));

    /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
    size_t nxout = grid->nx;
    size_t nyout = grid->ny;
    size_t nzout = grid->nz;
    float dxout = grid->dx;
    float dyout = grid->dy;
    float dzout = grid->dz;

    WRITE_HEADER_V0( dump_type::particle_dump, sp->id, sp->q/sp->m, fileIO, step, rank, nproc);

    dim[0] = sp->np;
    WRITE_ARRAY_HEADER( p_buf, 1, dim, fileIO );

    // Copy a PBUF_SIZE hunk of the particle list into the particle
    // buffer, timecenter it and write it out. This is done this way to
    // guarantee the particle list unchanged while not requiring too
    // much memory.

    // FIXME: WITH A PIPELINED CENTER_P, PBUF NOMINALLY SHOULD BE QUITE
    // LARGE.

    particle_t * sp_p = sp->p;      sp->p      = p_buf;
    int sp_np         = sp->np;     sp->np     = 0;
    int sp_max_np     = sp->max_np; sp->max_np = PBUF_SIZE;
    for( buf_start=0; buf_start<sp_np; buf_start += PBUF_SIZE ) {
        sp->np = sp_np-buf_start; if( sp->np > PBUF_SIZE ) sp->np = PBUF_SIZE;
        COPY( sp->p, &sp_p[buf_start], sp->np );
        center_p( sp, interpolator_array );
        fileIO.write( sp->p, sp->np );
    }
    sp->p      = sp_p;
    sp->np     = sp_np;
    sp->max_np = sp_max_np;

    if( fileIO.close() ) ERROR(("File close failed on dump particles!!!"));
}
void BinaryDump::dump_hydro(
        const char *fbase,
        int step,
        hydro_array_t* hydro_array,
        species_t* sp,
        interpolator_array_t* interpolator_array,
        grid_t* grid,
        int ftag
        )
{
    char fname[256];
    FileIO fileIO;
    int dim[3];

    if( !sp ) ERROR(( "Invalid species \"%s\"", sp->name ));

    clear_hydro_array( hydro_array );
    accumulate_hydro_p( hydro_array, sp, interpolator_array );
    synchronize_hydro_array( hydro_array );

    if( !fbase ) ERROR(( "Invalid filename" ));

    if( rank==0 )
        MESSAGE(("Dumping \"%s\" hydro fields to \"%s\"",sp->name,fbase));

    if( ftag ) sprintf( fname, "%s.%li.%i", fbase, (long)step, rank );
    else       sprintf( fname, "%s.%i", fbase, rank );
    FileIOStatus status = fileIO.open(fname, io_write);
    if( status==fail) ERROR(( "Could not open \"%s\".", fname ));

    /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
    size_t nxout = grid->nx;
    size_t nyout = grid->ny;
    size_t nzout = grid->nz;
    float dxout = grid->dx;
    float dyout = grid->dy;
    float dzout = grid->dz;

    WRITE_HEADER_V0(dump_type::hydro_dump, sp->id, sp->q/sp->m, fileIO, step, rank, nproc);

    dim[0] = grid->nx+2;
    dim[1] = grid->ny+2;
    dim[2] = grid->nz+2;
    WRITE_ARRAY_HEADER( hydro_array->h, 3, dim, fileIO );
    fileIO.write( hydro_array->h, dim[0]*dim[1]*dim[2] );
    if( fileIO.close() ) ERROR(( "File close failed on dump hydro!!!" ));
}

HDF5Dump::HDF5Dump(int _rank, int _nproc) : Dump_Strategy(_rank, _nproc) {
#ifdef HAS_FIELD_COMP
    field_type_id = H5Tcreate_fields();
#endif
#ifdef HAS_HYDRO_COMP
    hydro_type_id = H5Tcreate_hydro();
#endif
#ifdef HAS_PARTICLE_COMP
    particle_type_id = H5Tcreate_particle();
#endif
#ifdef HAS_EXPLICIT_ASYNC
    es_field = H5EScreate();
    es_hydro = H5EScreate();
    es_particle = H5EScreate();
#else
    es_field = es_hydro = es_particle = H5I_INVALID_HID;
#endif
    char *env_fprefix = getenv("VPIC_FILE_PREFIX");
    fprefix = (env_fprefix) ? env_fprefix : "";
}

HDF5Dump::~HDF5Dump() {
#ifdef HAS_EXPLICIT_ASYNC
    H5ESclose(es_field);
    H5ESclose(es_hydro);
    H5ESclose(es_particle);
#endif

#ifdef HAS_FIELD_COMP
    H5Tclose(field_type_id);
    #endif
#ifdef HAS_HYDRO_COMP
    H5Tclose(hydro_type_id);
#endif
#ifdef HAS_PARTICLE_COMP
    H5Tclose(particle_type_id);
#endif
}

hid_t HDF5Dump::H5Tcreate_fields(void) {
  hid_t comp_type = H5Tcreate(H5T_COMPOUND, sizeof(field_t));

  if (this->rank == 0)
      std::cout << "-- Field type size is " << sizeof(field_t) << std::endl;

  H5Tinsert(comp_type, "ex", HOFFSET(field_t, ex), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "ey", HOFFSET(field_t, ey), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "ez", HOFFSET(field_t, ez), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "div_e_err", HOFFSET(field_t, div_e_err),
            H5T_NATIVE_FLOAT);

  H5Tinsert(comp_type, "cbx", HOFFSET(field_t, cbx), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "cby", HOFFSET(field_t, cby), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "cbz", HOFFSET(field_t, cbz), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "div_b_err", HOFFSET(field_t, div_b_err),
            H5T_NATIVE_FLOAT);

  H5Tinsert(comp_type, "tcax", HOFFSET(field_t, tcax), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "tcay", HOFFSET(field_t, tcay), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "tcaz", HOFFSET(field_t, tcaz), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "rhob", HOFFSET(field_t, rhob), H5T_NATIVE_FLOAT);

  H5Tinsert(comp_type, "jfx", HOFFSET(field_t, jfx), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "jfy", HOFFSET(field_t, jfy), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "jfz", HOFFSET(field_t, jfz), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "rhof", HOFFSET(field_t, rhof), H5T_NATIVE_FLOAT);

  H5Tinsert(comp_type, "ematx", HOFFSET(field_t, ematx), H5T_NATIVE_SHORT);
  H5Tinsert(comp_type, "ematy", HOFFSET(field_t, ematy), H5T_NATIVE_SHORT);
  H5Tinsert(comp_type, "ematz", HOFFSET(field_t, ematz), H5T_NATIVE_SHORT);
  H5Tinsert(comp_type, "nmat", HOFFSET(field_t, nmat), H5T_NATIVE_SHORT);

  H5Tinsert(comp_type, "fmatx", HOFFSET(field_t, fmatx), H5T_NATIVE_SHORT);
  H5Tinsert(comp_type, "fmaty", HOFFSET(field_t, fmaty), H5T_NATIVE_SHORT);
  H5Tinsert(comp_type, "fmatz", HOFFSET(field_t, fmatz), H5T_NATIVE_SHORT);
  H5Tinsert(comp_type, "cmat", HOFFSET(field_t, cmat), H5T_NATIVE_SHORT);

  return comp_type;
}

hid_t HDF5Dump::H5Tcreate_hydro(void) {
  hid_t comp_type = H5Tcreate(H5T_COMPOUND, sizeof(hydro_t));

  if (this->rank == 0)
      std::cout << "-- Hydro type size is " << sizeof(hydro_t) << std::endl;

  H5Tinsert(comp_type, "jx", HOFFSET(hydro_t, jx), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "jy", HOFFSET(hydro_t, jy), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "jz", HOFFSET(hydro_t, jz), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "rho", HOFFSET(hydro_t, rho), H5T_NATIVE_FLOAT);

  H5Tinsert(comp_type, "px", HOFFSET(hydro_t, px), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "py", HOFFSET(hydro_t, py), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "pz", HOFFSET(hydro_t, pz), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "ke", HOFFSET(hydro_t, ke), H5T_NATIVE_FLOAT);

  H5Tinsert(comp_type, "txx", HOFFSET(hydro_t, txx), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "tyy", HOFFSET(hydro_t, tyy), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "tzz", HOFFSET(hydro_t, tzz), H5T_NATIVE_FLOAT);

  H5Tinsert(comp_type, "tyz", HOFFSET(hydro_t, tyz), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "tzx", HOFFSET(hydro_t, tzx), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "txy", HOFFSET(hydro_t, txy), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "pad", HOFFSET(hydro_t, _pad), H5T_NATIVE_DOUBLE);

  return comp_type;
}

hid_t HDF5Dump::H5Tcreate_particle(void) {
  hid_t comp_type = H5Tcreate(H5T_COMPOUND, sizeof(particle_t));

  if (this->rank == 0)
      std::cout << "-- Particle type size is " << sizeof(particle_t) << std::endl;

  H5Tinsert(comp_type, "dx", HOFFSET(particle_t, dx), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "dy", HOFFSET(particle_t, dy), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "dz", HOFFSET(particle_t, dz), H5T_NATIVE_FLOAT);

  H5Tinsert(comp_type, "i", HOFFSET(particle_t, i), H5T_NATIVE_INT);

  H5Tinsert(comp_type, "ux", HOFFSET(particle_t, ux), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "uy", HOFFSET(particle_t, uy), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "uz", HOFFSET(particle_t, uz), H5T_NATIVE_FLOAT);
  H5Tinsert(comp_type, "w", HOFFSET(particle_t, w), H5T_NATIVE_FLOAT);

  return comp_type;
}