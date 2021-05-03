#ifndef Dump_Strategy_h
#define Dump_Strategy_h

#include <iostream>
#include <unordered_map>
#include <vector>

//#define DUMP_INFO_DEBUG 1
//#define H5_ASYNC 1
#ifdef H5_ASYNC
#  include "h5_vol_external_async_native.h"
#endif // H5_ASYNC
// #define CHUNK_FLAG 1

//#define METADATA_COLL_WRITE 1
//#define TRUE 1

#define HAS_FIELD_COMP 1
#define HAS_PARTICLE_COMP 1
#define HAS_HYDRO_COMP 1

// #define HAS_FIELD_MAP 1
// #define HAS_PARTICLE_MAP 1
// #define HAS_HYDRO_MAP 1

//#define HAS_INDEPENDENT_IO 1

#include <cassert>
#include <mpi.h> // TODO: it would be good if this didn't have to know about MPI

// TODO: should I drop the ./src here?
#include "../field_advance/field_advance.h"
#include "../sf_interface/sf_interface.h"
#include "../species_advance/species_advance.h"
#include "../util/io/FileIO.h"
#include "../util/io/FileUtils.h"
#include "../util/util_base.h"

#include "dump.h"
#include "dumpmacros.h"

#ifdef VPIC_ENABLE_HDF5
#  include "hdf5.h"             // from the lib
#  include "hdf5_header_info.h" // from vpic
#  ifdef HAS_DAOS_VOL_EXT
#    include "daos_vol_public.h"
#  endif
#endif

#ifdef VPIC_ENABLE_OPENPMD
#  include <openPMD/openPMD.hpp>
#endif

//#define N_FILE_N_PROCESS 1
//#define TEST_MPIIO 1

// #define IO_LOG
#ifdef IO_LOG
#  define io_log(x)                                                            \
    do {                                                                       \
      if (rank == 0) {                                                         \
        std::cerr << __FILE__ << "(" << __func__ << ":" << __LINE__ << ")["    \
                  << rank << "]: " << x << std::endl;                          \
        std::cerr.flush();                                                     \
      }                                                                        \
    } while (0)
#else // IO_LOG
#  define io_log(x)                                                            \
    do {                                                                       \
    } while (0)
#endif // IO_LOG

#define fpp(x, y, z) f[VOXEL(x, y, z, grid->nx, grid->ny, grid->nz)]

#define DUMP_FIELD_TO_HDF5(DSET_NAME, ATTRIBUTE_NAME, ELEMENT_TYPE)            \
  {                                                                            \
    dset_id = H5Dcreate(group_id, DSET_NAME, ELEMENT_TYPE, filespace,          \
                        H5P_DEFAULT, dcpl_id, H5P_DEFAULT);                    \
    temp_buf_index = 0;                                                        \
    for (size_t i(1); i < grid->nx + 1; i++) {                                 \
      for (size_t j(1); j < grid->ny + 1; j++) {                               \
        for (size_t k(1); k < grid->nz + 1; k++) {                             \
          temp_buf[temp_buf_index] =                                           \
              FIELD_ARRAY_NAME->fpp(i, j, k).ATTRIBUTE_NAME;                   \
          temp_buf_index = temp_buf_index + 1;                                 \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    dataspace_id = H5Dget_space(dset_id);                                      \
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, global_offset, NULL,     \
                        global_count, NULL);                                   \
    H5Dwrite(dset_id, ELEMENT_TYPE, memspace, dataspace_id, plist_id,          \
             temp_buf);                                                        \
    H5Sclose(dataspace_id);                                                    \
    H5Dclose(dset_id);                                                         \
  }

// TODO: naming a macro so close to existing functions AND data is not a good
// define to do C-style indexing
#define _hydro(x, y, z)                                                        \
  hydro_array->h[VOXEL(x, y, z, grid->nx, grid->ny, grid->nz)]

#define WRITE_H5_FILE(group_id_p, data_buf_p, type_p, dname_p)                 \
  {                                                                            \
    hid_t dset_id = H5Dcreate(group_id_p, dname_p, type_p, filespace,          \
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);          \
    H5Dwrite(dset_id, type_p, memspace, filespace, io_plist_id, data_buf_p);   \
    H5Dclose(dset_id);                                                         \
  }

#define WRITE_MPI_FILE(dname_p, offset_p, data_buf_p, count_p, type_p)         \
  {                                                                            \
    MPI_File fh;                                                               \
    MPI_Status status;                                                         \
    sprintf(fname, "%s%s/%s_%d_%s.h5", fprefix, subparticle_scratch, sp->name, \
            step, dname_p);                                                    \
    if (mpi_rank == 0)                                                         \
      printf("fname= %s \n", fname);                                           \
    MPI_Info info;                                                             \
    MPI_Info_create(&info);                                                    \
    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,    \
                  info, &fh);                                                  \
    MPI_File_write_at(fh, offset_p, data_buf_p, count_p, type_p, &status);     \
    MPI_Info_free(&info);                                                      \
    MPI_File_close(&fh);                                                       \
  }

#define DUMP_HYDRO_TO_HDF5(DSET_NAME, ATTRIBUTE_NAME, ELEMENT_TYPE)            \
  {                                                                            \
    dset_id = H5Dcreate(group_id, DSET_NAME, ELEMENT_TYPE, filespace,          \
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);                \
    temp_buf_index = 0;                                                        \
    for (size_t i(1); i < grid->nx + 1; i++) {                                 \
      for (size_t j(1); j < grid->ny + 1; j++) {                               \
        for (size_t k(1); k < grid->nz + 1; k++) {                             \
          temp_buf[temp_buf_index] = _hydro(i, j, k).ATTRIBUTE_NAME;           \
          temp_buf_index = temp_buf_index + 1;                                 \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    dataspace_id = H5Dget_space(dset_id);                                      \
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, global_offset, NULL,     \
                        global_count, NULL);                                   \
    H5Dwrite(dset_id, ELEMENT_TYPE, memspace, dataspace_id, plist_id,          \
             temp_buf);                                                        \
    H5Sclose(dataspace_id);                                                    \
    H5Dclose(dset_id);                                                         \
  }

// Runtime inheritance is obviously not very "VPIC like", as we will [probably]
// incur a penalty for the vtable lookup, but given we're about to do IO this
// is very negligible.
class Dump_Strategy {
public:
  int rank, nproc, num_step;

  Dump_Strategy(int _rank, int _nproc) : rank(_rank), nproc(_nproc) {} // empty

  virtual ~Dump_Strategy(){};

  virtual void dump_fields(const char *fbase, int step, grid_t *grid,
                           field_array_t *field_array, int ftag) = 0;
  virtual void dump_hydro(const char *fbase, int step,
                          hydro_array_t *hydro_array, species_t *sp,
                          interpolator_array_t *interpolator_array,
                          grid_t *grid, int ftag) = 0;
  virtual void dump_particles(const char *fbase, species_t *sp, grid_t *grid,
                              int step,
                              interpolator_array_t *interpolator_array,
                              int ftag) = 0;
};

class BinaryDump : public Dump_Strategy {
public:
  using Dump_Strategy::Dump_Strategy; // inherit constructor
  // BinaryDump(int _rank, int _nproc ) : Dump_Strategy(_rank, _nproc ){ } //
  // empty

  // TODO: now we pass rank and step, ftag has odd semanticds
  void dump_fields(const char *fbase, int step, grid_t *grid,
                   field_array_t *field_array, int ftag);
  void dump_hydro(const char *fbase, int step, hydro_array_t *hydro_array,
                  species_t *sp, interpolator_array_t *interpolator_array,
                  grid_t *grid, int ftag);
  void dump_particles(const char *fbase, species_t *sp, grid_t *grid, int step,
                      interpolator_array_t *interpolator_array, int ftag);
};

#ifdef VPIC_ENABLE_HDF5

struct field_dump_flag_t {
  bool ex = true, ey = true, ez = true, div_e_err = true;
  bool cbx = true, cby = true, cbz = true, div_b_err = true;
  bool tcax = true, tcay = true, tcaz = true, rhob = true;
  bool jfx = true, jfy = true, jfz = true, rhof = true;
  bool ematx = true, ematy = true, ematz = true, nmat = true;
  bool fmatx = true, fmaty = true, fmatz = true, cmat = true;
  void disableE() { ex = false, ey = false, ez = false, div_e_err = false; }

  void disableCB() { cbx = false, cby = false, cbz = false, div_b_err = false; }

  void disableTCA() { tcax = false, tcay = false, tcaz = false, rhob = false; }

  void disableJF() { jfx = false, jfy = false, jfz = false, rhof = false; }

  void disableEMAT() {
    ematx = false, ematy = false, ematz = false, nmat = false;
  }

  void disableFMAT() {
    fmatx = false, fmaty = false, fmatz = false, cmat = false;
  }

  void resetToDefaults() {
    ex = true, ey = true, ez = true, div_e_err = true;
    cbx = true, cby = true, cbz = true, div_b_err = true;
    tcax = true, tcay = true, tcaz = true, rhob = true;
    jfx = true, jfy = true, jfz = true, rhof = true;
    ematx = true, ematy = true, ematz = true, nmat = true;
    fmatx = true, fmaty = true, fmatz = true, cmat = true;
  }

  bool enabledE() { return ex && ey && ez; }

  bool enabledCB() { return cbx && cby && cbz; }

  bool enabledTCA() { return tcax && tcay && tcaz; }

  bool enabledJF() { return jfx && jfy && jfz; }

  bool enabledEMAT() { return ematx && ematy && ematz; }

  bool enabledFMAT() { return fmatx && fmaty && fmatz; }
};

struct hydro_dump_flag_t {
  bool jx = true, jy = true, jz = true, rho = true;
  bool px = true, py = true, pz = true, ke = true;
  bool txx = true, tyy = true, tzz = true;
  bool tyz = true, tzx = true, txy = true;

  void disableJ() { jx = false, jy = false, jz = false, rho = false; }

  void disableP() { px = false, py = false, pz = false, ke = false; }

  void disableTD() // Stress diagonal
  {
    txx = false, tyy = false, tzz = false;
  }

  void disableTOD() // Stress off-diagonal
  {
    tyz = false, tzx = false, txy = false;
  }
  void resetToDefaults() {
    jx = true, jy = true, jz = true, rho = true;
    px = true, py = true, pz = true, ke = true;
    txx = true, tyy = true, tzz = true;
    tyz = true, tzx = true, txy = true;
  }

  bool enabledJ() { return jx && jy && jz; }

  bool enabledP() { return px && py && pz; }

  bool enabledTD() { return txx && tyy && tzz; }

  bool enabledTOD() { return tyz && tzx && txy; }
};

class HDF5Dump : public Dump_Strategy {
  std::unordered_map<species_id, size_t> tframe_map;

private:
#  ifdef HAS_FIELD_COMP
  field_t *field_buf;
  hid_t field_type_id;
#  endif
#  ifdef HAS_HYDRO_COMP
  hydro_t *hydro_buf;
  hid_t hydro_type_id;
#  endif
#  ifdef HAS_PARTICLE_COMP
  hid_t particle_type_id;
#  endif
  hid_t es_field, es_hydro, es_particle;
  const char *fprefix;
  const char *dump_dir;
  int indep_meta;
  int async;
  double io_time;

  hid_t H5Tcreate_fields(void);
  hid_t H5Tcreate_hydro(void);
  hid_t H5Tcreate_particle(void);

  inline hid_t H5Fcreate_wrap(const char *filename, unsigned flags,
                              hid_t fcpl_id, hid_t fapl_id, hid_t es_id) {
    hid_t file_id;
    double t = MPI_Wtime();
    if (async)
      file_id = H5Fcreate_async(filename, flags, fcpl_id, fapl_id, es_id);
    else
      file_id = H5Fcreate(filename, flags, fcpl_id, fapl_id);
    io_time += (MPI_Wtime() - t);
    return file_id;
  }

  inline herr_t H5Fclose_wrap(hid_t file_id, hid_t es_id) {
    herr_t err;
    double t = MPI_Wtime();
    if (async)
      err = H5Fclose_async(file_id, es_id);
    else
      err = H5Fclose(file_id);
    io_time += (MPI_Wtime() - t);
    return err;
  }

  inline hid_t H5Gcreate_wrap(hid_t loc_id, const char *name, hid_t lcpl_id,
                              hid_t gcpl_id, hid_t gapl_id, hid_t es_id) {
    hid_t grp_id;
    double t = MPI_Wtime();
    if (async)
      grp_id = H5Gcreate_async(loc_id, name, lcpl_id, gcpl_id, gapl_id, es_id);
    else
      grp_id = H5Gcreate2(loc_id, name, lcpl_id, gcpl_id, gapl_id);
    io_time += (MPI_Wtime() - t);
    return grp_id;
  }

  inline herr_t H5Gclose_wrap(hid_t group_id, hid_t es_id) {
    herr_t err;
    double t = MPI_Wtime();
    if (async)
      err = H5Gclose_async(group_id, es_id);
    else
      err = H5Gclose(group_id);
    io_time += (MPI_Wtime() - t);
    return err;
  }

  inline hid_t H5Dcreate_wrap(hid_t loc_id, const char *name, hid_t type_id,
                              hid_t space_id, hid_t lcpl_id, hid_t dcpl_id,
                              hid_t dapl_id, hid_t es_id) {
    hid_t dset_id;
    double t = MPI_Wtime();
    if (async)
      dset_id = H5Dcreate_async(loc_id, name, type_id, space_id, lcpl_id,
                                dcpl_id, dapl_id, es_id);
    else
      dset_id =
          H5Dcreate(loc_id, name, type_id, space_id, lcpl_id, dcpl_id, dapl_id);
    io_time += (MPI_Wtime() - t);
    return dset_id;
  }

  inline herr_t H5Dwrite_wrap(hid_t dset_id, hid_t mem_type_id,
                              hid_t mem_space_id, hid_t file_space_id,
                              hid_t dxpl_id, const void *buf, hid_t es_id) {
    herr_t err;
    double t = MPI_Wtime();
    if (async)
      err = H5Dwrite_async(dset_id, mem_type_id, mem_space_id, file_space_id,
                           dxpl_id, buf, es_id);
    else
      err = H5Dwrite(dset_id, mem_type_id, mem_space_id, file_space_id, dxpl_id,
                     buf);
    io_time += (MPI_Wtime() - t);
    return err;
  }

  inline herr_t H5Dclose_wrap(hid_t dset_id, hid_t es_id) {
    herr_t err;
    double t = MPI_Wtime();
    if (async)
      err = H5Dclose_async(dset_id, es_id);
    else
      err = H5Dclose(dset_id);
    io_time += (MPI_Wtime() - t);
    return err;
  }

  inline hid_t H5Mcreate_wrap(hid_t loc_id, const char *name, hid_t key_type_id,
                              hid_t val_type_id, hid_t lcpl_id, hid_t mcpl_id,
                              hid_t mapl_id, hid_t es_id) {
    hid_t map_id;
    double t = MPI_Wtime();
    if (async)
      map_id = H5Mcreate_async(loc_id, name, key_type_id, val_type_id, lcpl_id,
                               mcpl_id, mapl_id, es_id);
    else
      map_id = H5Mcreate(loc_id, name, key_type_id, val_type_id, lcpl_id,
                         mcpl_id, mapl_id);
    io_time += (MPI_Wtime() - t);
    return map_id;
  }

  inline herr_t H5Mput_wrap(hid_t map_id, hid_t key_mem_type_id,
                            const void *key, hid_t val_mem_type_id,
                            const void *value, hid_t dxpl_id, hid_t es_id) {
    herr_t err;
    double t = MPI_Wtime();
    if (async)
      err = H5Mput_async(map_id, key_mem_type_id, key, val_mem_type_id, value,
                         dxpl_id, es_id);
    else
      err =
          H5Mput(map_id, key_mem_type_id, key, val_mem_type_id, value, dxpl_id);
    io_time += (MPI_Wtime() - t);
    return err;
  }

  inline herr_t H5Mclose_wrap(hid_t map_id, hid_t es_id) {
    herr_t err;
    double t = MPI_Wtime();
    if (async)
      err = H5Mclose_async(map_id, es_id);
    else
      err = H5Mclose(map_id);
    io_time += (MPI_Wtime() - t);
    return err;
  }

  inline void asyncWait(hid_t es_id, uint64_t timeout) {
    size_t num_in_progress = 0;
    hbool_t err_occurred = 0;
    double t = MPI_Wtime();

    if (rank == 0)
      std::cout << "Here in wait" << std::endl;

    /* check if all operations in event set have completed */
    H5ESwait(es_id, timeout, &num_in_progress, &err_occurred);
    if (err_occurred || ((timeout == H5ES_WAIT_FOREVER) && num_in_progress)) {
      ERROR(("Failed to complete field async I/O \n"));
    }
    io_time += (MPI_Wtime() - t);
  }

public:
  HDF5Dump(int _rank, int _nproc);
  virtual ~HDF5Dump();

  // TODO: replace these with a common dump interface
  // Declare vars to use
  hydro_dump_flag_t hydro_dump_flag;
  field_dump_flag_t field_dump_flag;

  /**
   * @brief Dump field data to the HDf5 file
   *         Author: Bin Dong  dbin@lbl.gov
   *         https://crd.lbl.gov/bin-dong
   *         Nov 2020
   * @param fbase
   * @param step
   * @param grid
   * @param field_array
   * @param ftag
   */
  void dump_fields(const char *fbase, int step, grid_t *grid,
                   field_array_t *field_array, int ftag) {
    // double dump_field_uptime = uptime();
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    double t_start = uptime();

#  ifdef DUMP_INFO_DEBUG
    printf("MPI rank = %d, size = %d \n", mpi_rank, mpi_size);
    // printf("base dir for field: %s \n", fdParams.baseDir);
    // printf("stride x y z  = (%ld, %ld, %ld)\n", fdParams.stride_x,
    // fdParams.stride_y, fdParams.stride_z);
    printf("grid x, y z  = (%d, %d, %d) \n", grid->nx, grid->ny, grid->nz);
    printf("domain loc (x0, y0, z0) -> (x1, y1, z1) = (%f, %f, %f) -> (%f, %f, "
           "%f) \n",
           grid->x0, grid->y0, grid->z0, grid->x1, grid->y1, grid->z1);
    // printf("global->topology_x, y, z =  %f, %f, %f \n ", global->topology_x,
    // global->topology_y, global->topology_z);
    printf("grid -> sx, sy, sz =  (%d, %d, %d), nv=%d \n", grid->sx, grid->sy,
           grid->sz, grid->nv);
#  endif // DUMP_INFO_DEBUG

    char fname[1024];
    char field_scratch[256];
    char subfield_scratch[512];

    if (strcmp(dump_dir, ".") == 0)
      sprintf(field_scratch, "%s", "field_hdf5");
    else
      sprintf(field_scratch, "%s/%s", dump_dir, "field_hdf5");

    // FileUtils::makeDirectory(field_scratch);
    sprintf(subfield_scratch, "%s_T_%d", field_scratch, step);
    // FileUtils::makeDirectory(subfield_scratch);

    sprintf(fname, "%s%s_%s_%d.h5", fprefix, subfield_scratch, "fields", step);
    // double el1 = uptime();

    //    int file_exist(const char *filename)
    //{
    //    struct stat buffer;
    //    return (stat(filename, &buffer) == 0);
    //}

    // struct stat buffer;
    // if((stat(fname, &buffer) == 0)){
    //    file_exist_flag  = 1;
    //    if(!mpi_rank)
    //        printf("Write original files /w HDF5! \n");
    // }
    // file_exist_flag = 0;

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    // H5Pset_alignment(plist_id, 4194304, 4194304);
    /*if(H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST) <
      0){ exit(-1);
      }*/

#  ifdef METADATA_COLL_WRITE
    if (!mpi_rank)
      printf("Enable collective metadata write !\n");
    H5Pset_coll_metadata_write(plist_id, TRUE);
#  endif // METADATA_COLL_WRITE

#  ifdef HAS_DAOS_VOL_EXT
    if (indep_meta)
      H5daos_set_all_ind_metadata_ops(plist_id, 1);
#  endif

    hid_t file_id =
        H5Fcreate_wrap(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id, es_field);
    H5Pclose(plist_id);

#  ifdef HAS_DAOS_VOL_EXT
    if (indep_meta)
      sprintf(fname, "Timestep_%d_%d", step, rank);
    else
#  endif
      sprintf(fname, "Timestep_%d", step);

    hid_t group_id = H5Gcreate_wrap(file_id, fname, H5P_DEFAULT, H5P_DEFAULT,
                                    H5P_DEFAULT, es_field);

    // io_log("TimeHDF5Open:  " << uptime() - el1
    //  << " s"); // Easy to handle results for scripts
    // double el2 = uptime();

    /*
    // Create a variable list of field values to output.
    size_t numvars = std::min(global->fdParams.output_vars.bitsum(),
    total_field_variables); size_t * varlist = new size_t[numvars];

    for(size_t i(0), c(0); i<total_field_variables; i++)
    if(global->fdParams.output_vars.bitset(i)) varlist[c++] = i;

    printf("\nBEGIN_OUTPUT: numvars = %zd \n", numvars);*/

    /*
       typedef struct field {
       float ex,   ey,   ez,   div_e_err;     // Electric field and div E error
       float cbx,  cby,  cbz,  div_b_err;     // Magnetic field and div B error
       float tcax, tcay, tcaz, rhob;          // TCA fields and bound charge
       density float jfx,  jfy,  jfz,  rhof;          // Free current and charge
       density material_id ematx, ematy, ematz, nmat; // Material at edge
       centers and nodes material_id fmatx, fmaty, fmatz, cmat; // Material at
       face and cell centers } field_t;*/
    // Local voxel mesh resolution.  Voxels are
    // indexed FORTRAN style 0:nx+1,0:ny+1,0:nz+1
    // with voxels 1:nx,1:ny,1:nz being non-ghost
    // voxels.

    /*
       typedef struct field {
       float ex,   ey,   ez,   div_e_err;     // Electric field and div E
       error float cbx,  cby,  cbz,  div_b_err;     // Magnetic field and div
       B error float tcax, tcay, tcaz, rhob;          // TCA fields and bound
       charge density float jfx,  jfy,  jfz,  rhof;          // Free current
       and charge density material_id ematx, ematy, ematz, nmat; // Material
       at edge centers and nodes material_id fmatx, fmaty, fmatz, cmat; //
       Material at face and cell centers } field_t;*/

#  ifdef HAS_FIELD_MAP
    hid_t map_id =
        H5Mcreate_wrap(group_id, "field", H5T_NATIVE_INT, field_type_id,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_field);

    for (int i(1); i < grid->nx + 1; i++)
      for (int j(1); j < grid->ny + 1; j++)
        for (int k(1); k < grid->nz + 1; k++) {
          int global_index = VOXEL(i, j, k, grid->nx, grid->ny, grid->nz);
          H5Mput_wrap(map_id, H5T_NATIVE_INT, &global_index, field_type_id,
                      &FIELD_ARRAY_NAME->fpp(i, j, k), H5P_DEFAULT, es_field);
        }

    H5Mclose_wrap(map_id, es_field);
#  else // HAS_FIELD_MAP
    // char  *field_var_name[] =
    // {"ex","ey","ez","div_e_err","cbx","cby","cbz","div_b_err","tcax","tcay","tcaz","rhob","jfx","jfy","jfz","rhof"};
    // Comment out for test only

    plist_id = H5Pcreate(H5P_DATASET_XFER);
#    ifdef HAS_INDEPENDENT_IO
    if (!mpi_rank)
      printf("\n ###\n VPIC Independent I/O! \n ###\n");
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
#    else  // HAS_INDEPENDENT_IO
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#    endif // HAS_INDEPENDENT_IO

    // H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *) &offset, NULL,
    // (hsize_t *) &numparticles, NULL);

    // global->topology_x

    hsize_t field_global_size[3], field_local_size[3], global_offset[3],
        global_count[3];
    field_global_size[0] = (grid->nx * grid->gpx);
    field_global_size[1] = (grid->ny * grid->gpy);
    field_global_size[2] = (grid->nz * grid->gpz);

    field_local_size[0] = grid->nx;
    field_local_size[1] = grid->ny;
    field_local_size[2] = grid->nz;

    // Convert rank to local decomposition
    int rx, ry, rz;
    UNVOXEL(mpi_rank, rx, ry, rz, grid->gpx, grid->gpy, grid->gpz);

    int mpi_rank_x, mpi_rank_y, mpi_rank_z;
    mpi_rank_x = rx;
    mpi_rank_y = ry;
    mpi_rank_z = rz;

    global_offset[0] = (grid->nx) * mpi_rank_x;
    global_offset[1] = (grid->ny) * mpi_rank_y;
    global_offset[2] = (grid->nz) * mpi_rank_z;

    global_count[0] = (grid->nx);
    global_count[1] = (grid->ny);
    global_count[2] = (grid->nz);

#    ifdef DUMP_INFO_DEBUG
    if (mpi_rank < 4) {
      printf("grid nx, ny nz  = (%d, %d, %d) \n", grid->nx, grid->ny, grid->nz);
      printf("mpi-rank = %d, rank index = (%d, %d, %d) \n", mpi_rank,
             mpi_rank_x, mpi_rank_y, mpi_rank_z);
      printf("global size   = %llu  %llu %llu \n", field_global_size[0],
             field_global_size[1], field_global_size[2]);
      printf("global_offset = %llu %llu %llu \n", global_offset[0],
             global_offset[1], global_offset[2]);
      printf("global_count  = %llu  %llu %llu \n", global_count[0],
             global_count[1], global_count[2]);
      fflush(stdout);
    }
#    endif // DUMP_INFO_DEBUG

#    ifdef CHUNK_FLAG
    // hsize_t chunk_dims[3] = {288, 24, 24};
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk_dims[3] = {grid->nx, grid->ny, grid->nz};

    H5Pset_chunk(dcpl_id, 3, chunk_dims);
    if (!mpi_rank)
      printf("Enable chunking !\n");
#    else  // CHUNK_FLAG
    hid_t dcpl_id = H5P_DEFAULT;
#    endif // CHUNK_FLAG

    hsize_t temp_buf_index;
    hid_t dset_id;
    hid_t filespace = H5Screate_simple(3, field_global_size, NULL);
    hid_t memspace = H5Screate_simple(3, field_local_size, NULL);

#    ifdef HAS_FIELD_COMP
    temp_buf_index = 0;

    dset_id = H5Dcreate_wrap(group_id, "field", field_type_id, filespace,
                             H5P_DEFAULT, dcpl_id, H5P_DEFAULT, es_field);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, global_offset, NULL,
                        global_count, NULL);

    if (async) {
      if (rank == 0)
        std::cout << "About to wait" << std::endl;
      asyncWait(es_field, H5ES_WAIT_FOREVER);
      if (rank == 0)
        std::cout << "Exiting wait" << std::endl;
    }

    if (field_buf)
      free(field_buf);

    field_buf = (field_t *)malloc(sizeof(field_t) * (grid->nx) * (grid->ny) *
                                  (grid->nz));

    for (int i(1); i < grid->nx + 1; i++) {
      for (int j(1); j < grid->ny + 1; j++) {
        for (int k(1); k < grid->nz + 1; k++) {
          field_buf[temp_buf_index] = FIELD_ARRAY_NAME->fpp(i, j, k);
          temp_buf_index++;
        }
      }
    }

    H5Dwrite_wrap(dset_id, field_type_id, memspace, filespace, plist_id,
                  field_buf, es_field);
    H5Dclose_wrap(dset_id, es_field);
#    else  // HAS_FIELD_COMP
    float *temp_buf =
        (float *)malloc(sizeof(float) * (grid->nx) * (grid->ny) * (grid->nz));

    if (field_dump_flag.ex)
      DUMP_FIELD_TO_HDF5("ex", ex, H5T_NATIVE_FLOAT);
    if (field_dump_flag.ey)
      DUMP_FIELD_TO_HDF5("ey", ey, H5T_NATIVE_FLOAT);
    if (field_dump_flag.ez)
      DUMP_FIELD_TO_HDF5("ez", ez, H5T_NATIVE_FLOAT);
    if (field_dump_flag.div_e_err)
      DUMP_FIELD_TO_HDF5("div_e_err", div_e_err, H5T_NATIVE_FLOAT);

    if (field_dump_flag.cbx)
      DUMP_FIELD_TO_HDF5("cbx", cbx, H5T_NATIVE_FLOAT);
    if (field_dump_flag.cby)
      DUMP_FIELD_TO_HDF5("cby", cby, H5T_NATIVE_FLOAT);
    if (field_dump_flag.cbz)
      DUMP_FIELD_TO_HDF5("cbz", cbz, H5T_NATIVE_FLOAT);
    if (field_dump_flag.div_b_err)
      DUMP_FIELD_TO_HDF5("div_b_err", div_b_err, H5T_NATIVE_FLOAT);

    if (field_dump_flag.tcax)
      DUMP_FIELD_TO_HDF5("tcax", tcax, H5T_NATIVE_FLOAT);
    if (field_dump_flag.tcay)
      DUMP_FIELD_TO_HDF5("tcay", tcay, H5T_NATIVE_FLOAT);
    if (field_dump_flag.tcaz)
      DUMP_FIELD_TO_HDF5("tcaz", tcaz, H5T_NATIVE_FLOAT);
    if (field_dump_flag.rhob)
      DUMP_FIELD_TO_HDF5("rhob", rhob, H5T_NATIVE_FLOAT);

    if (field_dump_flag.jfx)
      DUMP_FIELD_TO_HDF5("jfx", jfx, H5T_NATIVE_FLOAT);
    if (field_dump_flag.jfy)
      DUMP_FIELD_TO_HDF5("jfy", jfy, H5T_NATIVE_FLOAT);
    if (field_dump_flag.jfz)
      DUMP_FIELD_TO_HDF5("jfz", jfz, H5T_NATIVE_FLOAT);
    if (field_dump_flag.rhof)
      DUMP_FIELD_TO_HDF5("rhof", rhof, H5T_NATIVE_FLOAT);

    // H5T_NATIVE_SHORT  for material_id (typedef int16_t material_id)
    if (field_dump_flag.ematx)
      DUMP_FIELD_TO_HDF5("ematx", ematx, H5T_NATIVE_SHORT);
    if (field_dump_flag.ematy)
      DUMP_FIELD_TO_HDF5("ematy", ematy, H5T_NATIVE_SHORT);
    if (field_dump_flag.ematz)
      DUMP_FIELD_TO_HDF5("ematz", ematz, H5T_NATIVE_SHORT);
    if (field_dump_flag.nmat)
      DUMP_FIELD_TO_HDF5("nmat", nmat, H5T_NATIVE_SHORT);

    if (field_dump_flag.fmatx)
      DUMP_FIELD_TO_HDF5("fmatx", fmatx, H5T_NATIVE_SHORT);
    if (field_dump_flag.fmaty)
      DUMP_FIELD_TO_HDF5("fmaty", fmaty, H5T_NATIVE_SHORT);
    if (field_dump_flag.fmatz)
      DUMP_FIELD_TO_HDF5("fmatz", fmatz, H5T_NATIVE_SHORT);
    if (field_dump_flag.cmat)
      DUMP_FIELD_TO_HDF5("cmat", cmat, H5T_NATIVE_SHORT);

    free(temp_buf);
#    endif // HAS_FIELD_COMP
#    ifdef CHUNK_FLAG
    H5Pclose(dcpl_id);
#    endif // CHUNK_FLAG
    H5Pclose(plist_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
#  endif   // HAS_FIELD_MAP

    /*
       H5D_mpio_actual_io_mode_t actual_io_mode;
       H5Pget_mpio_actual_io_mode(plist_id, &actual_io_mode);

       switch(actual_io_mode){
       case H5D_MPIO_NO_COLLECTIVE:
       io_log("H5Pget_mpio_actual_io_mode: H5D_MPIO_NO_COLLECTIVE: ");
       break;
       case H5D_MPIO_CHUNK_INDEPENDENT:
       io_log("H5Pget_mpio_actual_io_mode: H5D_MPIO_CHUNK_INDEPENDENT: ");
       break;
       case H5D_MPIO_CHUNK_COLLECTIVE:
       io_log("H5Pget_mpio_actual_io_mode: H5D_MPIO_CHUNK_COLLECTIVE: ");
       break;
       case H5D_MPIO_CHUNK_MIXED:
       io_log("H5Pget_mpio_actual_io_mode: H5D_MPIO_CHUNK_MIXED: ");
       break;
       case H5D_MPIO_CONTIGUOUS_COLLECTIVE:
       io_log("H5Pget_mpio_actual_io_mode: H5D_MPIO_CONTIGUOUS_COLLECTIVE: ");
       break;
       default :
       io_log("H5Pget_mpio_actual_io_mode: None returend: ");
       break;
       }

       H5D_mpio_actual_chunk_opt_mode_t actual_chunk_opt_mode;
       H5Pget_mpio_actual_chunk_opt_mode(plist_id, &actual_chunk_opt_mode);
       switch(actual_chunk_opt_mode){
       case H5D_MPIO_NO_CHUNK_OPTIMIZATION:
       io_log("H5Pget_mpio_actual_chunk_opt_mode:
H5D_MPIO_NO_CHUNK_OPTIMIZATION: "); break; case H5D_MPIO_MULTI_CHUNK:
       io_log("H5Pget_mpio_actual_chunk_opt_mode: H5D_MPIO_MULTI_CHUNK: ");
       break;
    //  case H5D_MPIO_MULTI_CHUNK_NO_OPT:
    //      io_log("H5Pget_mpio_actual_chunk_opt_mode:
H5D_MPIO_MULTI_CHUNK_NO_OPT: ");
    //     break;
    case H5D_MPIO_LINK_CHUNK:
    io_log("H5Pget_mpio_actual_chunk_opt_mode: H5D_MPIO_LINK_CHUNK: ");
    break;
    default :
    io_log("H5Pget_mpio_actual_chunk_opt_mode: None returend: ");
    break;
    }

    uint32_t local_no_collective_cause,  global_no_collective_cause;
    H5Pget_mpio_no_collective_cause(plist_id, &local_no_collective_cause,
&global_no_collective_cause);

    switch(local_no_collective_cause){
    case H5D_MPIO_COLLECTIVE:
    io_log("local_no_collective_cause: H5D_MPIO_COLLECTIVE: ");
    break;
    case H5D_MPIO_SET_INDEPENDENT:
    io_log("local_no_collective_cause: H5D_MPIO_SET_INDEPENDENT: ");
    break;
    case H5D_MPIO_DATA_TRANSFORMS:
    io_log("local_no_collective_cause: H5D_MPIO_DATA_TRANSFORMS: ");
    break;
    //case H5D_MPIO_SET_MPIPOSIX:
    //    io_log("local_no_collective_cause: H5D_MPIO_SET_MPIPOSIX: ");
    //    break;
    case H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES:
    io_log("local_no_collective_cause: H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES:
"); break;
    //case H5D_MPIO_POINT_SELECTIONS:
    //    io_log("local_no_collective_cause: H5D_MPIO_POINT_SELECTIONS: ");
    //    break;
    // case H5D_MPIO_FILTERS:
    //    io_log("local_no_collective_cause: H5D_MPIO_FILTERS: ");
    //    break;
    default :
    io_log("local_no_collective_cause: None returend: ");
    break;
}


switch(global_no_collective_cause){
    case H5D_MPIO_COLLECTIVE:
        io_log("global_no_collective_cause: H5D_MPIO_COLLECTIVE: ");
        break;
    case H5D_MPIO_SET_INDEPENDENT:
        io_log("global_no_collective_cause: H5D_MPIO_SET_INDEPENDENT: ");
        break;
    case H5D_MPIO_DATA_TRANSFORMS:
        io_log("global_no_collective_cause: H5D_MPIO_DATA_TRANSFORMS: ");
        break;
        //case H5D_MPIO_SET_MPIPOSIX:
        //    io_log("global_no_collective_cause: H5D_MPIO_SET_MPIPOSIX: ");
        //    break;
    case H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES:
        io_log("global_no_collective_cause:
H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES: "); break;
        //case H5D_MPIO_POINT_SELECTIONS:
        //    io_log("global_no_collective_cause: H5D_MPIO_POINT_SELECTIONS: ");
        //    break;
        // case H5D_MPIO_FILTERS:
        //   io_log("global_no_collective_cause: H5D_MPIO_FILTERS: ");
        //   break;
    default :
        io_log("global_no_collective_cause: None returend: ");
        break;
}
*/

    // io_log("TimeHDF5Write: " << uptime() - el2 << " s");
    // double el3 = uptime();

    /*
    //Write metadata (geo original and geo dx/dy/dz) for ArrayUDF
    float attr_data[2][3];
    attr_data[0][0] = grid->x0;
    attr_data[0][1] = grid->y0;
    attr_data[0][2] = grid->z0;
    attr_data[1][0] = grid->dx;
    attr_data[1][1] = grid->dy;
    attr_data[1][2] = grid->dz;
    hsize_t dims[2];
    dims[0] = 2;
    dims[1] = 3;
    if(!file_exist_flag){
    hid_t va_geo_dataspace_id = H5Screate_simple(2, dims, NULL);
    hid_t va_geo_attribute_id = H5Acreate2(file_id, "VPIC-ArrayUDF-GEO",
    H5T_IEEE_F32BE, va_geo_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(va_geo_attribute_id, H5T_NATIVE_FLOAT, attr_data);
    H5Sclose(va_geo_dataspace_id);
    H5Aclose(va_geo_attribute_id);
    }
    */
    // if(!file_exist_flag)
    H5Gclose_wrap(group_id, es_field);
    H5Fclose_wrap(file_id, es_field);

    // io_log("TimeHDF5Close: " << uptime() - el3 << " s");

#  ifdef OUTPUT_XDMF
    if (mpi_rank == 0) {
      char const *output_xml_file = "./field_hdf5/hdf5_field.xdmf";
      char dimensions_3d[128];
      sprintf(dimensions_3d, "%llu %llu %llu", field_global_size[0],
              field_global_size[1], field_global_size[2]);
      char dimensions_4d[128];
      sprintf(dimensions_4d, "%llu %llu %llu %d", field_global_size[0],
              field_global_size[1], field_global_size[2], 3);
      char orignal[128];
      sprintf(orignal, "%f %f %f", grid->x0, grid->y0, grid->z0);
      char dxdydz[128];
      sprintf(dxdydz, "%f %f %f", grid->dx, grid->dy, grid->dz);

      // TODO: remove or let the user set
      int field_interval = 1;

      // TODO: remove this dependence on number of steps
      // std::cout << "num_step " << num_step << std::endl;

      int nframes = num_step / field_interval + 1;
      static int field_tframe = 0;

#    ifdef DUMP_INFO_DEBUG
      printf("         meta file : %s \n", output_xml_file);
      printf(" array dims per var: %s \n", dimensions_3d);
      printf("array dims all vars: %s \n", dimensions_4d);
      printf("            orignal: %s \n", orignal);
      printf("             dxdydz: %s \n", dxdydz);
      printf("            nframes: %d \n", nframes);
      printf("    field_interval: %d \n", field_interval);
      printf("       current step: %d \n", step);
      printf("       current step: %d \n", step);

      // printf("    Simulation time: %f \n", grid->t0);
      printf("             tframe: %d \n", field_tframe);
#    endif // DUMP_INFO_DEBUG

      // TODO: this footer dumping is more likely better done in a
      // destructor, rather than hoping a multiple division works out
      if (field_tframe >= 1) {
        if (field_tframe == (nframes - 1)) {
          invert_field_xml_item(output_xml_file, "fields", step, dimensions_4d,
                                dimensions_3d, 1);
        } else {
          invert_field_xml_item(output_xml_file, "fields", step, dimensions_4d,
                                dimensions_3d, 0);
        }
      } else {
        create_file_with_header(output_xml_file, dimensions_3d, orignal, dxdydz,
                                nframes, field_interval);
        if (field_tframe == (nframes - 1)) {
          invert_field_xml_item(output_xml_file, "fields", step, dimensions_4d,
                                dimensions_3d, 1);
        } else {
          invert_field_xml_item(output_xml_file, "fields", step, dimensions_4d,
                                dimensions_3d, 0);
        }
      }
      field_tframe++;
    }
#  endif // OUTPUT_XDMF
    double t_end = uptime();
    if (!rank)
      printf("Total dump field time for %d fields: %lf\n",
             (grid->nx) * (grid->ny) * (grid->nz), t_end - t_start);
  }
  /**
   * @brief dump_particles to the HDF5 file
   *         Author: Bin Dong  dbin@lbl.gov
   *            https://crd.lbl.gov/bin-dong
   *         Nov 2020
   * @param fbase
   * @param sp
   * @param grid
   * @param step
   * @param interpolator_array
   * @param ftag
   */
  void dump_particles(const char *fbase, species_t *sp, grid_t *grid, int step,
                      interpolator_array_t *interpolator_array, int ftag) {
    static int file_index = 0;
    file_index++;
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // double dump_particles_uptime = uptime();
    // time_t seconds = time(NULL);
    // printf("Atrank = %d, file_index = %d, dump_particles_uptime = %f,
    // epoch_seconds = %ld  \n ", mpi_rank, file_index, dump_particles_uptime,
    // seconds);

    char fname[1024];
    char group_name[1024];
    char particle_scratch[256];
    char subparticle_scratch[512];

    int np_local = sp->np;

    // get the total number of particles. in this example, output only electrons
    // sp = species_list;
    if (strcmp(dump_dir, ".") == 0)
      sprintf(particle_scratch, "%s", "particle_hdf5");
    else
      sprintf(particle_scratch, "%s/%s", dump_dir, "particle_hdf5");
    // FileUtils::makeDirectory(particle_scratch);
    sprintf(subparticle_scratch, "%s_T_%d", particle_scratch, step);
    // FileUtils::makeDirectory(subparticle_scratch);

    // TODO: Allow the user to set this
    int stride_particle_dump = 1;
    np_local = (sp->np + stride_particle_dump - 1) / stride_particle_dump;

    // make a copy of the part of particle data to be dumped
    // double ec1 = uptime();

    // particle_t *ALIGNED(128) p_buf = NULL;
    // if (!p_buf)
    //   MALLOC_ALIGNED(p_buf, np_local, 128);
    // particle_t *sp_p = sp->p;
    // sp->p = p_buf;
    // sp->np = np_local;
    // sp->max_np = np_local;

    // for (long long iptl = 0, i = 0; iptl < sp->np;
    //      iptl += stride_particle_dump, ++i) {
    //   COPY(&sp->p[i], &sp_p[iptl], 1);
    // }

    // center_p(sp, interpolator_array);

    // ec1 = uptime() - ec1;

#  ifndef N_FILE_N_PROCESS
    int np_local_max, np_local_min;
    MPI_Reduce(&np_local, &np_local_max, 1, MPI_INT, MPI_MAX, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&np_local, &np_local_min, 1, MPI_INT, MPI_MIN, 0,
               MPI_COMM_WORLD);
#  endif // N_FILE_N_PROCESS

    long long offset, numparticles = np_local;

#  ifndef N_FILE_N_PROCESS
    MPI_Scan(&numparticles, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    offset -= numparticles;
#  else  // N_FILE_N_PROCESS
    offset = 0;
#  endif // N_FILE_N_PROCESS

    // open HDF5 file in "particle/T.<step>/" subdirectory
    // filename: eparticle.h5p
#  ifndef N_FILE_N_PROCESS
    sprintf(fname, "%s%s_%s_%d.h5", fprefix, subparticle_scratch, sp->name,
            step);
#  else  // N_FILE_N_PROCESS
    sprintf(fname, "%s%s/%s_%d_p%d.h5", fprefix, subparticle_scratch, sp->name,
            step, mpi_rank);
#  endif // N_FILE_N_PROCESS

    // double el1 = uptime();

    // double t_start = uptime();

    hid_t file_plist_id = H5Pcreate(H5P_FILE_ACCESS);

#  ifndef N_FILE_N_PROCESS
    H5Pset_fapl_mpio(file_plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
#  endif // N_FILE_N_PROCESS
    H5Pset_mpi_params(file_plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

#  ifdef H5_ASYNC
    if (!mpi_rank)
      printf("Enable async on particle data");

    assert(H5Pset_vol_async(file_plist_id));
#  endif // H5_ASYNC

#  ifdef HAS_DAOS_VOL_EXT
    if (indep_meta) {
      H5daos_set_all_ind_metadata_ops(file_plist_id, 1);
    }
#  endif

    hid_t file_id = H5Fcreate_wrap(fname, H5F_ACC_TRUNC, H5P_DEFAULT,
                                   file_plist_id, es_particle);
    H5Pclose(file_plist_id);

#  ifdef HAS_DAOS_VOL_EXT
    if (indep_meta)
      sprintf(group_name, "/Timestep_%d_%d", step, rank);
    else
#  endif
      sprintf(group_name, "/Timestep_%d", step);

    double t_group_create1 = MPI_Wtime();

    hid_t group_id = H5Gcreate_wrap(file_id, group_name, H5P_DEFAULT,
                                    H5P_DEFAULT, H5P_DEFAULT, es_particle);

    double t_group_create2 = MPI_Wtime();
    if (!rank)
      printf("Time for group create operation: %lf\n",
             t_group_create2 - t_group_create1);

    // io_log("TimeHDF5Open:  " << uptime() - el1
    //  << " s"); // Easy to handle results for scripts
    // double el2 = uptime();

    {
      hid_t *part_group_ids = NULL;
      size_t n_groups = numparticles;
      long long total_particles;
      hid_t gcpl_id;

      if (!indep_meta) {
        // MPI_Allreduce(&numparticles, &n_groups, 1, MPI_LONG_LONG, MPI_SUM,
        //               MPI_COMM_WORLD);
        MPI_Bcast(&n_groups, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
      }

      gcpl_id = H5Pcreate(H5P_GROUP_CREATE);
      H5daos_set_object_class(gcpl_id, "SX");
      part_group_ids = (hid_t *) malloc(n_groups * sizeof(hid_t));

      MPI_Barrier(MPI_COMM_WORLD);
      if (rank == 0)
        printf("Creating %zu groups\n", n_groups);
      double t_group_create_part1 = MPI_Wtime();

      if (indep_meta) {
        for (long long i = 0; i < n_groups; i++) {
          char part_group_name[64];
          sprintf(part_group_name, "Part_%d", offset + i);

          part_group_ids[i] =
              H5Gcreate_wrap(group_id, part_group_name, H5P_DEFAULT,
                             gcpl_id, H5P_DEFAULT, es_particle);
        }
      } else {
        for (long long i = 0; i < n_groups; i++) {
          char part_group_name[64];
          sprintf(part_group_name, "Part_%d", i);

          if (rank == 0 && (i % 10000) == 0)
              printf(".");

          part_group_ids[i] =
              H5Gcreate_wrap(group_id, part_group_name, H5P_DEFAULT,
                             gcpl_id, H5P_DEFAULT, es_particle);
        }
        printf("\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);
      H5Pclose(gcpl_id);

      double t_group_create_part2 = MPI_Wtime();
      if (!rank)
        printf("Time for part group create operation: %lf\n",
               t_group_create_part2 - t_group_create_part1);

      for (long long i = 0; i < n_groups; i++)
        H5Gclose_wrap(part_group_ids[i], es_particle);
      free(part_group_ids);
    }
    MPI_Barrier(MPI_COMM_WORLD);

#  ifdef HAS_PARTICLE_MAP
    hid_t map_id =
        H5Mcreate_wrap(group_id, "particle", H5T_NATIVE_LLONG, particle_type_id,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_particle);

    if (!rank)
      printf("Writing %lld particles\n", sp->np);
    for (long long i(0); i < sp->np; i++) {
      long long particle_index = offset + i;
      H5Mput_wrap(map_id, H5T_NATIVE_LLONG, &particle_index, particle_type_id,
                  &sp->p[i], H5P_DEFAULT, es_particle);
    }

    H5Mclose_wrap(map_id, es_particle);
#  else // HAS_PARTICLE_MAP
    long long total_particles;

#    ifndef N_FILE_N_PROCESS
    MPI_Allreduce(&numparticles, &total_particles, 1, MPI_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
#    else  // N_FILE_N_PROCESS
    total_particles = np_local;
#    endif // N_FILE_N_PROCESS

    hid_t filespace = H5Screate_simple(1, (hsize_t *)&total_particles, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *)&offset, NULL,
                        (hsize_t *)&numparticles, NULL);

    // plist_id = H5P_DEFAULT;
    hid_t io_plist_id = H5Pcreate(H5P_DATASET_XFER);

#    ifndef N_FILE_N_PROCESS
#      ifdef HAS_INDEPENDENT_IO
    if (!mpi_rank) {
      printf("\n ###\n VPIC Independent I/O! \n ###\n");
    }
    H5Pset_dxpl_mpio(io_plist_id, H5FD_MPIO_INDEPENDENT);
#      else  // HAS_INDEPENDENT_IO
    // H5Pset_dxpl_mpio(io_plist_id, H5FD_MPIO_COLLECTIVE);
#      endif // HAS_INDEPENDENT_IO
#    endif   // N_FILE_N_PROCESS

#    ifdef H5_ASYNC
    H5Pset_dxpl_async(io_plist_id, true);
#    endif // H5_ASYNC
    hsize_t memspace_count_temp;
    hid_t memspace;
#    ifdef HAS_PARTICLE_COMP
    memspace_count_temp = numparticles;
    memspace = H5Screate_simple(1, &memspace_count_temp, NULL);
#    else  // HAS_PARTICLE_COMP
    memspace_count_temp = numparticles * 8;
    memspace = H5Screate_simple(1, &memspace_count_temp, NULL);
    hsize_t memspace_start = 0, memspace_stride = 8, memspace_count = np_local;
    // H5Sselect_hyperslab(memspace, H5S_SELECT_SET, &memspace_start,
    //                     &memspace_stride, &memspace_count, NULL);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, &memspace_start, NULL,
                        &memspace_count, NULL);
#    endif // HAS_PARTICLE_COMP

    // MPI_Info_set(info, "romio_cb_write", "disable");

#    ifdef HAS_PARTICLE_COMP
    hid_t dset_id =
        H5Dcreate_wrap(group_id, "particle", particle_type_id, filespace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_particle);

    double t_start = uptime();

    H5Dwrite_wrap(dset_id, particle_type_id, memspace, filespace, io_plist_id,
                  sp->p, es_particle);
    H5Dclose_wrap(dset_id, es_particle);
#    else // HAS_PARTICLE_COMP
    float *Pf = (float *)sp->p;
    int *Pi = (int *)sp->p;

#      ifdef TEST_MPIIO
    // Here we don't use the stripe but just for performance test
    if (!mpi_rank)
      printf("Test MPI-IO\n");
    WRITE_MPI_FILE("dX", offset * sizeof(float), Pf, numparticles, MPI_FLOAT);
    WRITE_MPI_FILE("dY", offset * sizeof(float), Pf, numparticles, MPI_FLOAT);
    WRITE_MPI_FILE("dZ", offset * sizeof(float), Pf, numparticles, MPI_FLOAT);
    WRITE_MPI_FILE("i", offset * sizeof(int), Pf, numparticles, MPI_INT);
    WRITE_MPI_FILE("ux", offset * sizeof(float), Pf, numparticles, MPI_FLOAT);
    WRITE_MPI_FILE("uy", offset * sizeof(float), Pf, numparticles, MPI_FLOAT);
    WRITE_MPI_FILE("uz", offset * sizeof(float), Pf, numparticles, MPI_FLOAT);
    WRITE_MPI_FILE("q", offset * sizeof(float), Pf, numparticles, MPI_FLOAT);
#      else // TEST_MPIIO
#        ifndef N_FILE_N_PROCESS
    if (!mpi_rank)
      printf("Test HDF5-IO Single \n");
#        else  // N_FILE_N_PROCESS
    if (!mpi_rank)
      printf("Test HDF5-IO N Files N Process\n");
#        endif // N_FILE_N_PROCESS
    // if(!mpi_rank )
    // io_log("++Particle Starting to write ) ");

    hid_t dx_id =
        H5Dcreate_wrap(group_id, "dX", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT, es_particle);
    hid_t dy_id =
        H5Dcreate_wrap(group_id, "dY", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT, es_particle);
    hid_t dz_id =
        H5Dcreate_wrap(group_id, "dZ", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT, es_particle);
    hid_t i_id =
        H5Dcreate_wrap(group_id, "i", H5T_NATIVE_INT, filespace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT, es_particle);
    hid_t ux_id =
        H5Dcreate_wrap(group_id, "ux", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT, es_particle);
    hid_t uy_id =
        H5Dcreate_wrap(group_id, "uy", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT, es_particle);
    hid_t uz_id =
        H5Dcreate_wrap(group_id, "uz", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT, es_particle);
    hid_t q_id =
        H5Dcreate_wrap(group_id, "q", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
                       H5P_DEFAULT, H5P_DEFAULT, es_particle);

    if (async) {
      asyncWait(es_particle, H5ES_WAIT_FOREVER);
    }

    double t_start = uptime();

    H5Dwrite_wrap(dx_id, H5T_NATIVE_FLOAT, memspace, filespace, io_plist_id, Pf,
                  es_particle);
    H5Dwrite_wrap(dy_id, H5T_NATIVE_FLOAT, memspace, filespace, io_plist_id,
                  Pf + 1, es_particle);
    H5Dwrite_wrap(dz_id, H5T_NATIVE_FLOAT, memspace, filespace, io_plist_id,
                  Pf + 2, es_particle);
    H5Dwrite_wrap(i_id, H5T_NATIVE_FLOAT, memspace, filespace, io_plist_id,
                  Pf + 3, es_particle);
    H5Dwrite_wrap(ux_id, H5T_NATIVE_FLOAT, memspace, filespace, io_plist_id,
                  Pf + 4, es_particle);
    H5Dwrite_wrap(uy_id, H5T_NATIVE_FLOAT, memspace, filespace, io_plist_id,
                  Pf + 5, es_particle);
    H5Dwrite_wrap(uz_id, H5T_NATIVE_FLOAT, memspace, filespace, io_plist_id,
                  Pf + 6, es_particle);
    H5Dwrite_wrap(q_id, H5T_NATIVE_FLOAT, memspace, filespace, io_plist_id,
                  Pf + 7, es_particle);

    H5Dclose_wrap(dx_id, es_particle);
    H5Dclose_wrap(dy_id, es_particle);
    H5Dclose_wrap(dz_id, es_particle);
    H5Dclose_wrap(i_id, es_particle);
    H5Dclose_wrap(ux_id, es_particle);
    H5Dclose_wrap(uy_id, es_particle);
    H5Dclose_wrap(uz_id, es_particle);
    H5Dclose_wrap(q_id, es_particle);

#      endif // TEST_MPIIO
#    endif   // HAS_PARTICLE_COMP
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Pclose(io_plist_id);
#  endif     // HAS_PARTICLE_MAP

    // io_log("TimeHDF5Write: " << uptime() - el2 << " s");
    // double el3 = uptime();

    H5Gclose_wrap(group_id, es_particle);
    H5Fclose_wrap(file_id, es_particle);

#  ifdef H5_ASYNC
    H5VLasync_finalize();
#  endif // H5_ASYNC

    if (async) {
      asyncWait(es_particle, H5ES_WAIT_FOREVER);
    }

    // io_log("TimeHDF5Close: " << uptime() - el3 << " s");
    double t_end = uptime();

    if (!rank)
      printf("(%d) Total dump %s particles time for %lld particles: %lf\n",
             step, sp->name, sp->np, t_end - t_start);
  }

  /**
   * @brief Dump hydro data to the HDf5 file
   *         Author: Bin Dong  dbin@lbl.gov
   *           https://crd.lbl.gov/bin-dong
   *         Nov 2020
   * @param fbase
   * @param step
   * @param hydro_array
   * @param sp
   * @param interpolator_array
   * @param grid
   * @param ftag
   */
  void dump_hydro(const char *fbase, int step, hydro_array_t *hydro_array,
                  species_t *sp, interpolator_array_t *interpolator_array,
                  grid_t *grid, int ftag) {
    // double dump_hydro_uptime = uptime();
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // double t_start = uptime();

    if (!sp) {
      ERROR(("Invalid species"));
    }

    clear_hydro_array(hydro_array);
    accumulate_hydro_p(hydro_array, sp, interpolator_array);
    synchronize_hydro_array(hydro_array);

    char hname[1024];
    char hydro_scratch[256];
    char subhydro_scratch[512];

    if (strcmp(dump_dir, ".") == 0)
      sprintf(hydro_scratch, "%s", "hydro_hdf5");
    else
      sprintf(hydro_scratch, "%s/%s", dump_dir, "hydro_hdf5");
    // FileUtils::makeDirectory(hydro_scratch);
    sprintf(subhydro_scratch, "%s_T_%d", hydro_scratch, step);
    // FileUtils::makeDirectory(subhydro_scratch);

    sprintf(hname, "%s%s_hydro_%s_%d.h5", fprefix, subhydro_scratch, sp->name,
            step);
    // double el1 = uptime();
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    /*
       if(H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST) <
       0){ exit(-1);
       }*/
    // if((fid = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id)) < 0)
    //    ERROR_RETURN;

    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

#  ifdef HAS_DAOS_VOL_EXT
    if (indep_meta)
      H5daos_set_all_ind_metadata_ops(plist_id, 1);
#  endif

    hid_t file_id =
        H5Fcreate_wrap(hname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id, es_hydro);
    H5Pclose(plist_id);

#  ifdef HAS_DAOS_VOL_EXT
    if (indep_meta)
      sprintf(hname, "Timestep_%d_%d", step, rank);
    else
#  endif
      sprintf(hname, "Timestep_%d", step);

    hid_t group_id = H5Gcreate_wrap(file_id, hname, H5P_DEFAULT, H5P_DEFAULT,
                                    H5P_DEFAULT, es_hydro);

    // el1 = uptime() - el1;
    // io_log("TimeHDF5Open:  " << el1
    //  << " s"); // Easy to handle results for scripts
    // double el2 = uptime();

    // Create a variable list of field values to output.
    // size_t numvars = std::min(global->fdParams.output_vars.bitsum(),
    // total_field_variables); size_t *varlist = new size_t[numvars];

    // for (size_t i(0), c(0); i < total_field_variables; i++)
    //    if (global->fdParams.output_vars.bitset(i))
    //        varlist[c++] = i;

    // printf("\nBEGIN_OUTPUT: numvars = %zu \n", numvars);

    // typedef struct hydro {
    //  float jx, jy, jz, rho; // Current and charge density => <q v_i f>, <q f>
    //  float px, py, pz, ke;  // Momentum and K.E. density  => <p_i f>, <m c^2
    //  (gamma-1) f> float txx, tyy, tzz;   // Stress diagonal            =>
    //  <p_i v_j f>, i==j float tyz, tzx, txy;   // Stress off-diagonal => <p_i
    //  v_j f>, i!=j float _pad[2];         // 16-byte align
    //} hydro_t;

    // typedef struct hydro_array {
    //  hydro_t * ALIGNED(128) h;
    //  grid_t * g;
    //} hydro_array_t;

#  ifdef HAS_FIELD_MAP
    hid_t map_id =
        H5Mcreate_wrap(group_id, "hydro", H5T_NATIVE_INT, hydro_type_id,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_field);

    for (int i(1); i < grid->nx + 1; i++)
      for (int j(1); j < grid->ny + 1; j++)
        for (int k(1); k < grid->nz + 1; k++) {
          int global_index = VOXEL(i, j, k, grid->nx, grid->ny, grid->nz);
          H5Mput_wrap(map_id, H5T_NATIVE_INT, &global_index, field_type_id,
                      &_hydro(i, j, k), H5P_DEFAULT, es_field);
        }

    H5Mclose_wrap(map_id, es_field);
#  else
    // char  *field_var_name[] =
    // {"ex","ey","ez","div_e_err","cbx","cby","cbz","div_b_err","tcax","tcay","tcaz","rhob","jfx","jfy","jfz","rhof"};
    // Comment out for test only

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *) &offset, NULL,
    // (hsize_t *) &numparticles, NULL);

    // global->topology_x

    hsize_t hydro_global_size[3], hydro_local_size[3], global_offset[3],
        global_count[3];
    hydro_global_size[0] = (grid->nx * grid->gpx);
    hydro_global_size[1] = (grid->ny * grid->gpy);
    hydro_global_size[2] = (grid->nz * grid->gpz);

    hydro_local_size[0] = grid->nx;
    hydro_local_size[1] = grid->ny;
    hydro_local_size[2] = grid->nz;

    int mpi_rank_x, mpi_rank_y, mpi_rank_z;
    UNVOXEL(mpi_rank, mpi_rank_x, mpi_rank_y, mpi_rank_z, grid->gpx, grid->gpy,
            grid->gpz);

    global_offset[0] = (grid->nx) * mpi_rank_x;
    global_offset[1] = (grid->ny) * mpi_rank_y;
    global_offset[2] = (grid->nz) * mpi_rank_z;

    global_count[0] = (grid->nx);
    global_count[1] = (grid->ny);
    global_count[2] = (grid->nz);

#    ifdef DUMP_INFO_DEBUG
    printf("global size   = %llu %llu %llu \n", hydro_global_size[0],
           hydro_global_size[1], hydro_global_size[2]);
    printf("global_offset = %llu %llu %llu \n", global_offset[0],
           global_offset[1], global_offset[2]);
    printf("global_count  = %llu %llu %llu \n", global_count[0],
           global_count[1], global_count[2]);
    printf("mpi-rank = %d, rank index = (%d, %d, %d) \n", mpi_rank, mpi_rank_x,
           mpi_rank_y, mpi_rank_z);
    fflush(stdout);
#    endif

    hid_t filespace = H5Screate_simple(3, hydro_global_size, NULL);
    hid_t memspace = H5Screate_simple(3, hydro_local_size, NULL);

    // typedef struct hydro {
    //  float jx, jy, jz, rho; // Current and charge density => <q v_i f>, <q f>
    //  float px, py, pz, ke;  // Momentum and K.E. density  => <p_i f>, <m c^2
    //  (gamma-1) f> float txx, tyy, tzz;   // Stress diagonal            =>
    //  <p_i v_j f>, i==j float tyz, tzx, txy;   // Stress off-diagonal => <p_i
    //  v_j f>, i!=j float _pad[2];         // 16-byte align
    //} hydro_t;

    hsize_t temp_buf_index;
    hid_t dset_id;
#    ifdef HAS_HYDRO_COMP
    dset_id = H5Dcreate_wrap(group_id, "hydro", hydro_type_id, filespace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_hydro);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, global_offset, NULL,
                        global_count, NULL);

    if (async) {
      asyncWait(es_hydro, H5ES_WAIT_FOREVER);
    }

    if (hydro_buf)
      free(hydro_buf);

    hydro_buf = (hydro_t *)malloc(sizeof(hydro_t) * (grid->nx) * (grid->ny) *
                                  (grid->nz));

    temp_buf_index = 0;
    for (int i(1); i < grid->nx + 1; i++) {
      for (int j(1); j < grid->ny + 1; j++) {
        for (int k(1); k < grid->nz + 1; k++) {
          hydro_buf[temp_buf_index] = _hydro(i, j, k);
          temp_buf_index = temp_buf_index + 1;
        }
      }
    }

    H5Dwrite_wrap(dset_id, hydro_type_id, memspace, filespace, plist_id,
                  hydro_buf, es_hydro);
    H5Dclose_wrap(dset_id, es_hydro);
#    else
    float *temp_buf =
        (float *)malloc(sizeof(float) * (grid->nx) * (grid->ny) * (grid->nz));

    if (hydro_dump_flag.jx)
      DUMP_HYDRO_TO_HDF5("jx", jx, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.jy)
      DUMP_HYDRO_TO_HDF5("jy", jy, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.jz)
      DUMP_HYDRO_TO_HDF5("jz", jz, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.rho)
      DUMP_HYDRO_TO_HDF5("rho", rho, H5T_NATIVE_FLOAT);

    if (hydro_dump_flag.px)
      DUMP_HYDRO_TO_HDF5("px", px, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.py)
      DUMP_HYDRO_TO_HDF5("py", py, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.pz)
      DUMP_HYDRO_TO_HDF5("pz", pz, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.ke)
      DUMP_HYDRO_TO_HDF5("ke", ke, H5T_NATIVE_FLOAT);

    if (hydro_dump_flag.txx)
      DUMP_HYDRO_TO_HDF5("txx", txx, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.tyy)
      DUMP_HYDRO_TO_HDF5("tyy", tyy, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.tzz)
      DUMP_HYDRO_TO_HDF5("tzz", tzz, H5T_NATIVE_FLOAT);

    if (hydro_dump_flag.tyz)
      DUMP_HYDRO_TO_HDF5("tyz", tyz, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.tzx)
      DUMP_HYDRO_TO_HDF5("tzx", tzx, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.txy)
      DUMP_HYDRO_TO_HDF5("txy", txy, H5T_NATIVE_FLOAT);

    free(temp_buf);
#    endif
    H5Pclose(plist_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
#  endif

    // el2 = uptime() - el2;
    // io_log("TimeHDF5Write: " << el2 << " s");

    // double el3 = uptime();

    // Write metadata (geo original and geo dx/dy/dz) for ArrayUDF
    /*
       float attr_data[2][3];
       attr_data[0][0] = grid->x0;
       attr_data[0][1] = grid->y0;
       attr_data[0][2] = grid->z0;
       attr_data[1][0] = grid->dx;
       attr_data[1][1] = grid->dy;
       attr_data[1][2] = grid->dz;
       hsize_t dims[2];
       dims[0] = 2;
       dims[1] = 3;
       hid_t va_geo_dataspace_id = H5Screate_simple(2, dims, NULL);
       hid_t va_geo_attribute_id = H5Acreate2(file_id, "VPIC-ArrayUDF-GEO",
       H5T_IEEE_F32BE, va_geo_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
       H5Awrite(va_geo_attribute_id, H5T_NATIVE_FLOAT, attr_data);
       H5Sclose(va_geo_dataspace_id);
       H5Aclose(va_geo_attribute_id);*/

    H5Gclose_wrap(group_id, es_hydro);
    H5Fclose_wrap(file_id, es_hydro);

    // el3 = uptime() - el3;
    // io_log("TimeHDF5Close: " << el3 << " s");

#  ifdef OUTPUT_XDMF
    if (mpi_rank == 0) {
      char output_xml_file[128];
      sprintf(output_xml_file, "./%s/%s%s%s", "hydro_hdf5", "hydro-", sp->name,
              ".xdmf");
      char dimensions_3d[128];
      sprintf(dimensions_3d, "%llu %llu %llu", hydro_global_size[0],
              hydro_global_size[1], hydro_global_size[2]);
      char dimensions_4d[128];
      sprintf(dimensions_4d, "%llu %llu %llu %d", hydro_global_size[0],
              hydro_global_size[1], hydro_global_size[2], 3);
      char orignal[128];
      sprintf(orignal, "%f %f %f", grid->x0, grid->y0, grid->z0);
      char dxdydz[128];
      sprintf(dxdydz, "%f %f %f", grid->dx, grid->dy, grid->dz);

      // TODO: remove or let user set
      int hydro_interval = 1;

      // TODO: remove this dependence on number of steps
      int nframes = num_step / hydro_interval + 1;

      const int tframe = tframe_map[sp->id];

#    ifdef DUMP_INFO_DEBUG
      printf("         meta file : %s \n", output_xml_file);
      printf(" array dims per var: %s \n", dimensions_3d);
      printf("array dims all vars: %s \n", dimensions_4d);
      printf("            orignal: %s \n", orignal);
      printf("             dxdydz: %s \n", dxdydz);
      printf("            nframes: %d \n", nframes);
      printf("    hydro_fields_interval: %d \n", hydro_interval);
      printf("       current step: %d \n", step);
      printf("    Simulation time: %f \n", grid->t0);
      printf("             tframe: %d \n", tframe);
#    endif // DUMP_INFO_DEBUG

      // TODO: why doesnt this just use the cstr?
      char speciesname_new[128];
      sprintf(speciesname_new, "hydro_%s", sp->name);
      if (tframe >= 1) {
        if (tframe == (nframes - 1)) {
          invert_hydro_xml_item(output_xml_file, speciesname_new, step,
                                dimensions_4d, dimensions_3d, 1);
        } else {
          invert_hydro_xml_item(output_xml_file, speciesname_new, step,
                                dimensions_4d, dimensions_3d, 0);
        }
      } else {
        create_file_with_header(output_xml_file, dimensions_3d, orignal, dxdydz,
                                nframes, hydro_interval);
        if (tframe == (nframes - 1)) {
          invert_hydro_xml_item(output_xml_file, speciesname_new, step,
                                dimensions_4d, dimensions_3d, 1);
        } else {
          invert_hydro_xml_item(output_xml_file, speciesname_new, step,
                                dimensions_4d, dimensions_3d, 0);
        }
      }
      tframe_map[sp->id]++;
    }
#  endif // OUTPUT_XDMF
    //   double t_end = uptime();
    //   if (!rank)
    //     printf("Total dump %s hydro time for %d fields: %lf\n", sp->name,
    //            (grid->nx) * (grid->ny) * (grid->nz), t_end - t_start);
  }
};
#endif // VPIC_ENABLE_HDF5

#ifdef VPIC_ENABLE_OPENPMD
class OpenPMDDump : public Dump_Strategy {
public:
  // openPMD::Series* series;
  using Dump_Strategy::Dump_Strategy; // inherit constructor

  // std::string file_type = ".h5";
  std::string file_type = ".bp";

  void dump_fields(const char *fbase, int step, grid_t *grid,
                   field_array_t *field_array, int ftag) {
    std::cout << "Writing openPMD data" << std::endl;

    std::string full_file_name = fbase + file_type;

    // if (series == nullptr) {
    std::cout << "init series" << std::endl;
    openPMD::Series series = openPMD::Series(
        full_file_name, openPMD::AccessType::CREATE, MPI_COMM_WORLD);
    //}

    std::cout << "Writing iteration " << step << std::endl;
    auto i = series.iterations[step];
    // TODO: it would be nice to set these...
    // series.setAuthor( "Axel Huebl <a.huebl@hzdr.de>");
    // series.setMachine( "Hall Probe 5000, Model 3");
    i.setAttribute("vacuum", true);

    auto cB = i.meshes["B"];
    auto E = i.meshes["E"];
    auto J = i.meshes["J"];
    auto Tca = i.meshes["Tca"];
    auto Emat = i.meshes["Emat"];
    auto Fmat = i.meshes["Fmat"];
    auto Rho = i.meshes["Rho"];
    auto DivErr = i.meshes["DivErr"];

    // record components
    auto Cbx = cB["x"];
    auto Cby = cB["y"];
    auto Cbz = cB["z"];

    auto Ex = E["x"];
    auto Ey = E["y"];
    auto Ez = E["z"];

    auto Jx = J["x"];
    auto Jy = J["y"];
    auto Jz = J["z"];

    auto Tcax = Tca["x"];
    auto Tcay = Tca["y"];
    auto Tcaz = Tca["z"];

    auto Ematx = Emat["x"];
    auto Ematy = Emat["y"];
    auto Ematz = Emat["z"];

    auto Fmatx = Fmat["x"];
    auto Fmaty = Fmat["y"];
    auto Fmatz = Fmat["z"];

    auto RhoB = Rho["B"];
    auto RhoF = Rho["F"];

    auto DivEErr = DivErr["E"];
    auto DivBErr = DivErr["B"];

    // TODO: set unitDimension so the anaylsis software knows what fields
    // things are
    //
    // // TODO: add timers for the convert and for the write

    size_t gnx = (grid->nx * grid->gpx);
    size_t gny = (grid->ny * grid->gpy);
    size_t gnz = (grid->nz * grid->gpz);
    openPMD::Extent global_extent = {gnx, gny, gnz};

    openPMD::Datatype datatype = openPMD::determineDatatype<float>();
    openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

    Cbx.resetDataset(dataset);
    Cby.resetDataset(dataset);
    Cbz.resetDataset(dataset);

    Ex.resetDataset(dataset);
    Ey.resetDataset(dataset);
    Ez.resetDataset(dataset);

    Jx.resetDataset(dataset);
    Jy.resetDataset(dataset);
    Jz.resetDataset(dataset);

    Tcax.resetDataset(dataset);
    Tcay.resetDataset(dataset);
    Tcaz.resetDataset(dataset);

    Ematx.resetDataset(dataset);
    Ematy.resetDataset(dataset);
    Ematz.resetDataset(dataset);

    Fmatx.resetDataset(dataset);
    Fmaty.resetDataset(dataset);
    Fmatz.resetDataset(dataset);

    RhoB.resetDataset(dataset);
    RhoF.resetDataset(dataset);

    DivEErr.resetDataset(dataset);
    DivBErr.resetDataset(dataset);

    // TODO: hoist this conversion code, as is it used elsewhere
    // Convert rank to local x/y/z
    int rx, ry, rz;
    UNVOXEL(rank, rx, ry, rz, grid->gpx, grid->gpy, grid->gpz);

    size_t nx = grid->nx;
    size_t ny = grid->ny;
    size_t nz = grid->nz;

    // NOTE: this assumes a static mesh decomposition in nx/ny/nz
    size_t global_offset_x = (nx)*rx;
    size_t global_offset_y = (ny)*ry;
    size_t global_offset_z = (nz)*rz;

    openPMD::Offset chunk_offset = {global_offset_x, global_offset_y,
                                    global_offset_z};
    openPMD::Extent chunk_extent = {nx, ny, nz};

    std::cout << "Local offset "
              << " x: " << global_offset_x << " y: " << global_offset_y
              << " z: " << global_offset_z << std::endl;

    // Store a local copy of the data which we pull out of the AoS
    std::vector<float> cbx_data;
    std::vector<float> cby_data;
    std::vector<float> cbz_data;

    std::vector<float> ex_data;
    std::vector<float> ey_data;
    std::vector<float> ez_data;

    std::vector<float> jx_data;
    std::vector<float> jy_data;
    std::vector<float> jz_data;

    std::vector<float> tcax_data;
    std::vector<float> tcay_data;
    std::vector<float> tcaz_data;

    // TODO: these are material_id (ints not floats)
    std::vector<float> ematx_data;
    std::vector<float> ematy_data;
    std::vector<float> ematz_data;

    std::vector<float> fmatx_data;
    std::vector<float> fmaty_data;
    std::vector<float> fmatz_data;
    // end todo

    std::vector<float> rhob_data;
    std::vector<float> rhof_data;

    std::vector<float> divb_data;
    std::vector<float> dive_data;

    size_t nv = nx * ny * nz;

    // TODO: resize here will zero out the data which we don't need, we
    // could change to a different semantic to avoid this
    cbx_data.resize(nv);
    cby_data.resize(nv);
    cbz_data.resize(nv);

    ex_data.resize(nv);
    ey_data.resize(nv);
    ez_data.resize(nv);

    jx_data.resize(nv);
    jy_data.resize(nv);
    jz_data.resize(nv);

    tcax_data.resize(nv);
    tcay_data.resize(nv);
    tcaz_data.resize(nv);

    ematx_data.resize(nv);
    ematy_data.resize(nv);
    ematz_data.resize(nv);

    fmatx_data.resize(nv);
    fmaty_data.resize(nv);
    fmatz_data.resize(nv);

    rhob_data.resize(nv);
    rhof_data.resize(nv);

    divb_data.resize(nv);
    dive_data.resize(nv);

    // TODO: make this AoS to SoA conversion a function

    // We could do 1D here, but we don't really care about the ghosts, and
    // we can thread over nz/ny (collapsed?) Go over non-ghosts and grab
    // just that data into a dense array
    for (size_t k = 1; k < grid->nz + 1; k++) {
      for (size_t j = 1; j < grid->ny + 1; j++) {
        for (size_t i = 1; i < grid->nx + 1; i++) {
          int local_index = VOXEL(i - 1, j - 1, k - 1, grid->nx - 2,
                                  grid->ny - 2, grid->nz - 2);
          int global_index = VOXEL(i, j, k, grid->nx, grid->ny, grid->nz);

          cbx_data[local_index] = field_array->f[global_index].cbx;
          cby_data[local_index] = field_array->f[global_index].cby;
          cbz_data[local_index] = field_array->f[global_index].cbz;

          ex_data[local_index] = field_array->f[global_index].ex;
          ey_data[local_index] = field_array->f[global_index].ey;
          ez_data[local_index] = field_array->f[global_index].ez;

          jx_data[local_index] = field_array->f[global_index].jfx;
          jy_data[local_index] = field_array->f[global_index].jfy;
          jz_data[local_index] = field_array->f[global_index].jfz;

          tcax_data[local_index] = field_array->f[global_index].tcax;
          tcay_data[local_index] = field_array->f[global_index].tcay;
          tcaz_data[local_index] = field_array->f[global_index].tcaz;

          ematx_data[local_index] = field_array->f[global_index].ematx;
          ematy_data[local_index] = field_array->f[global_index].ematy;
          ematz_data[local_index] = field_array->f[global_index].ematz;

          fmatx_data[local_index] = field_array->f[global_index].fmatx;
          fmaty_data[local_index] = field_array->f[global_index].fmaty;
          fmatz_data[local_index] = field_array->f[global_index].fmatz;

          rhob_data[local_index] = field_array->f[global_index].rhob;
          rhof_data[local_index] = field_array->f[global_index].rhof;

          dive_data[local_index] = field_array->f[global_index].div_e_err;
          divb_data[local_index] = field_array->f[global_index].div_b_err;
        }
      }
    }

    Cbx.storeChunk(cbx_data, chunk_offset, chunk_extent);
    Cby.storeChunk(cby_data, chunk_offset, chunk_extent);
    Cbz.storeChunk(cbz_data, chunk_offset, chunk_extent);

    Ex.storeChunk(ex_data, chunk_offset, chunk_extent);
    Ey.storeChunk(ey_data, chunk_offset, chunk_extent);
    Ez.storeChunk(ez_data, chunk_offset, chunk_extent);

    Jx.storeChunk(jx_data, chunk_offset, chunk_extent);
    Jy.storeChunk(jy_data, chunk_offset, chunk_extent);
    Jz.storeChunk(jz_data, chunk_offset, chunk_extent);

    Tcax.storeChunk(tcax_data, chunk_offset, chunk_extent);
    Tcay.storeChunk(tcay_data, chunk_offset, chunk_extent);
    Tcaz.storeChunk(tcaz_data, chunk_offset, chunk_extent);

    Ematx.storeChunk(ematx_data, chunk_offset, chunk_extent);
    Ematy.storeChunk(ematy_data, chunk_offset, chunk_extent);
    Ematz.storeChunk(ematz_data, chunk_offset, chunk_extent);

    Fmatx.storeChunk(fmatx_data, chunk_offset, chunk_extent);
    Fmaty.storeChunk(fmaty_data, chunk_offset, chunk_extent);
    Fmatz.storeChunk(fmatz_data, chunk_offset, chunk_extent);

    RhoB.storeChunk(rhob_data, chunk_offset, chunk_extent);
    RhoF.storeChunk(rhof_data, chunk_offset, chunk_extent);

    DivEErr.storeChunk(dive_data, chunk_offset, chunk_extent);
    DivBErr.storeChunk(divb_data, chunk_offset, chunk_extent);

    series.flush();
  }

  void dump_particles(const char *fbase, species_t *sp, grid_t *grid, int step,
                      interpolator_array_t *interpolator_array, int ftag) {
    std::string full_file_name = fbase + file_type;

    std::cout << "writing particles to " << full_file_name << std::endl;

    // if (series == nullptr) {
    openPMD::Series series = openPMD::Series(
        full_file_name, openPMD::AccessType::CREATE, MPI_COMM_WORLD);
    //}

    auto i = series.iterations[step];

    // TODO: set these
    i.setTime((float)step);
    i.setDt(1.0);
    i.setTimeUnitSI(1.0);

    auto &p = i.particles[sp->name];

    const int np = sp->np;

    // TODO: this could be a function call as it's used elsewhere (in hdf5)
    unsigned long long total_particles, offset;
    unsigned long long numparticles = np;
    MPI_Allreduce(&numparticles, &total_particles, 1, MPI_LONG_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Scan(&numparticles, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    offset -= numparticles;

    openPMD::Extent global_extent = {total_particles};
    openPMD::Datatype datatype = openPMD::determineDatatype<float>();
    openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

    auto px = p["position"]["x"];
    auto pxo = p["positionOffset"]["x"];

    auto py = p["position"]["y"];
    auto pyo = p["positionOffset"]["y"];

    auto pz = p["position"]["z"];
    auto pzo = p["positionOffset"]["z"];

    auto ux = p["velocity"]["x"];
    auto uy = p["velocity"]["y"];
    auto uz = p["velocity"]["z"];

    px.resetDataset(dataset);
    pxo.resetDataset(dataset);

    py.resetDataset(dataset);
    pyo.resetDataset(dataset);

    pz.resetDataset(dataset);
    pzo.resetDataset(dataset);

    ux.resetDataset(dataset);
    uy.resetDataset(dataset);
    uz.resetDataset(dataset);
    // convert data to SoA, allowing the user to chunk the operation

    // TODO: Add code the convert to global offsets
#  ifndef PMD_MAX_IO_CHUNK             // in particles
#    define PMD_MAX_IO_CHUNK 16777216; // 512MB total write
#  endif
    const int max_chunk = PMD_MAX_IO_CHUNK;

    // Loop over all particles in chunks
    for (int i = 0; i < np; i += max_chunk) {
      // We have to be careful as the last chunk may not be full
      // Find how many are left and do that many
      size_t to_write = std::min(np - i, max_chunk);

      // Convert the chunk ready to write
      std::vector<float> x_pos;
      std::vector<float> x_off;
      x_pos.resize(to_write);
      x_off.resize(to_write);

      std::vector<float> y_pos;
      std::vector<float> y_off;
      y_pos.resize(to_write);
      y_off.resize(to_write);

      std::vector<float> z_pos;
      std::vector<float> z_off;
      z_pos.resize(to_write);
      z_off.resize(to_write);

      std::vector<float> ux_pos;
      ux_pos.resize(to_write);

      std::vector<float> uy_pos;
      uy_pos.resize(to_write);

      std::vector<float> uz_pos;
      uz_pos.resize(to_write);

      for (int j = 0; j < to_write; j++) {
        // TODO: do I need to center the particles?
        auto &particle = sp->p[i + j];

        x_pos[j] = particle.dx;
        y_pos[j] = particle.dy;
        z_pos[j] = particle.dz;

        ux_pos[j] = particle.ux;
        uy_pos[j] = particle.uy;
        uz_pos[j] = particle.uz;

        std::array<int, 4> gi = global_particle_index(particle.i, grid, rank);
        x_off[j] = (float)gi[1];
        y_off[j] = (float)gi[2];
        z_off[j] = (float)gi[3];
      }

      // Base offset plus i to account for chunks
      auto o = openPMD::Offset{offset + i};
      auto e = openPMD::Extent{to_write};
      px.storeChunk(x_pos, o, e);
      pxo.storeChunk(x_off, o, e);

      py.storeChunk(y_pos, o, e);
      pyo.storeChunk(y_off, o, e);

      pz.storeChunk(z_pos, o, e);
      pzo.storeChunk(z_off, o, e);

      ux.storeChunk(ux_pos, o, e);
      uy.storeChunk(uy_pos, o, e);
      uz.storeChunk(uz_pos, o, e);

      series.flush();
    }
  }
  void dump_hydro(const char *fbase, int step, hydro_array_t *hydro_array,
                  species_t *sp, interpolator_array_t *interpolator_array,
                  grid_t *grid, int ftag) {
    std::string full_file_name = fbase + file_type;

    std::cout << "OpenPMD dumping hydro to " << full_file_name << std::endl;

    // if (series == nullptr) {
    openPMD::Series series = openPMD::Series(
        full_file_name, openPMD::AccessType::CREATE, MPI_COMM_WORLD);
    //}

    auto i = series.iterations[step];

    // TODO: set these
    i.setTime((float)step);
    i.setDt(1.0);
    i.setTimeUnitSI(1.0);

    if (!sp)
      ERROR(("Invalid species \"%s\"", sp->name));

    // TODO: do we want each backend to have to explicitly call these
    // manually? Or, as it is common, should we hoist it to the VPIC
    // call-site
    clear_hydro_array(hydro_array);
    accumulate_hydro_p(hydro_array, sp, interpolator_array);
    synchronize_hydro_array(hydro_array);

    if (!fbase)
      ERROR(("Invalid filename"));

    if (rank == 0)
      MESSAGE(("Dumping \"%s\" hydro fields to \"%s\"", sp->name, fbase));

    // Write data
    // float jx, jy, jz, rho; // Current and charge density => <q v_i f>, <q
    // f> float px, py, pz, ke;  // Momentum and K.E. density  => <p_i f>,
    // <m c^2 (gamma-1) f> float txx, tyy, tzz;   // Stress diagonal => <p_i
    // v_j f>, i==j float tyz, tzx, txy;   // Stress off-diagonal        =>
    // <p_i v_j f>, i!=j
    auto J = i.meshes["J"];
    auto P = i.meshes["P"];
    auto T = i.meshes["T"];
    auto _Ke = i.meshes["Ke"];
    auto _Rho = i.meshes["Rho"];

    auto Jx = J["x"];
    auto Jy = J["y"];
    auto Jz = J["z"];

    auto Px = P["x"];
    auto Py = P["y"];
    auto Pz = P["z"];

    auto Txx = T["xx"];
    auto Tyy = T["yy"];
    auto Tzz = T["zz"];
    auto Tyz = T["yz"];
    auto Tzx = T["zx"];
    auto Txy = T["xy"];

    auto Rho = _Rho["rho"]; // TODO: bad name..
    auto Ke = _Ke["ke"];    // TODO: bad name..

    size_t gnx = (grid->nx * grid->gpx);
    size_t gny = (grid->ny * grid->gpy);
    size_t gnz = (grid->nz * grid->gpz);
    openPMD::Extent global_extent = {gnx, gny, gnz};

    openPMD::Datatype datatype = openPMD::determineDatatype<float>();
    openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

    Jx.resetDataset(dataset);
    Jy.resetDataset(dataset);
    Jz.resetDataset(dataset);

    Px.resetDataset(dataset);
    Py.resetDataset(dataset);
    Pz.resetDataset(dataset);

    Txx.resetDataset(dataset);
    Tyy.resetDataset(dataset);
    Tzz.resetDataset(dataset);
    Tyz.resetDataset(dataset);
    Tzx.resetDataset(dataset);
    Txy.resetDataset(dataset);

    Rho.resetDataset(dataset);
    Ke.resetDataset(dataset);

    // TODO: hoist this conversion code, as is it used elsewhere
    // Convert rank to local x/y/z
    int rx, ry, rz;
    UNVOXEL(rank, rx, ry, rz, grid->gpx, grid->gpy, grid->gpz);

    size_t nx = grid->nx;
    size_t ny = grid->ny;
    size_t nz = grid->nz;

    // NOTE: this assumes a static mesh decomposition in nx/ny/nz
    size_t global_offset_x = (nx)*rx;
    size_t global_offset_y = (ny)*ry;
    size_t global_offset_z = (nz)*rz;

    openPMD::Offset chunk_offset = {global_offset_x, global_offset_y,
                                    global_offset_z};
    openPMD::Extent chunk_extent = {nx, ny, nz};

    std::cout << "Local offset "
              << " x: " << global_offset_x << " y: " << global_offset_y
              << " z: " << global_offset_z << std::endl;

    std::vector<float> jx_data;
    std::vector<float> jy_data;
    std::vector<float> jz_data;

    std::vector<float> px_data;
    std::vector<float> py_data;
    std::vector<float> pz_data;

    std::vector<float> txx_data;
    std::vector<float> tyy_data;
    std::vector<float> tzz_data;
    std::vector<float> tyz_data;
    std::vector<float> tzx_data;
    std::vector<float> txy_data;

    std::vector<float> rho_data;
    std::vector<float> ke_data;

    size_t nv = nx * ny * nz;

    jx_data.resize(nv);
    jy_data.resize(nv);
    jz_data.resize(nv);

    px_data.resize(nv);
    py_data.resize(nv);
    pz_data.resize(nv);

    txx_data.resize(nv);
    tyy_data.resize(nv);
    tzz_data.resize(nv);
    tyz_data.resize(nv);
    tzx_data.resize(nv);
    txy_data.resize(nv);

    rho_data.resize(nv);
    ke_data.resize(nv);

    // Transpose AoS to SoAs
    for (size_t k = 1; k < grid->nz + 1; k++) {
      for (size_t j = 1; j < grid->ny + 1; j++) {
        for (size_t i = 1; i < grid->nx + 1; i++) {
          int local_index = VOXEL(i - 1, j - 1, k - 1, grid->nx - 2,
                                  grid->ny - 2, grid->nz - 2);
          int global_index = VOXEL(i, j, k, grid->nx, grid->ny, grid->nz);

          jx_data[local_index] = hydro_array->h[global_index].jx;
          jy_data[local_index] = hydro_array->h[global_index].jy;
          jz_data[local_index] = hydro_array->h[global_index].jz;

          px_data[local_index] = hydro_array->h[global_index].px;
          py_data[local_index] = hydro_array->h[global_index].py;
          pz_data[local_index] = hydro_array->h[global_index].pz;

          txx_data[local_index] = hydro_array->h[global_index].txx;
          tyy_data[local_index] = hydro_array->h[global_index].tyy;
          tzz_data[local_index] = hydro_array->h[global_index].tzz;
          tyz_data[local_index] = hydro_array->h[global_index].tyz;
          tzx_data[local_index] = hydro_array->h[global_index].tzx;
          txy_data[local_index] = hydro_array->h[global_index].txy;

          rho_data[local_index] = hydro_array->h[global_index].rho;
          ke_data[local_index] = hydro_array->h[global_index].ke;
        }
      }
    }

    Jx.storeChunk(jx_data, chunk_offset, chunk_extent);
    Jy.storeChunk(jy_data, chunk_offset, chunk_extent);
    Jz.storeChunk(jz_data, chunk_offset, chunk_extent);

    Px.storeChunk(px_data, chunk_offset, chunk_extent);
    Py.storeChunk(py_data, chunk_offset, chunk_extent);
    Pz.storeChunk(pz_data, chunk_offset, chunk_extent);

    Txx.storeChunk(txx_data, chunk_offset, chunk_extent);
    Tyy.storeChunk(tyy_data, chunk_offset, chunk_extent);
    Tzz.storeChunk(tzz_data, chunk_offset, chunk_extent);
    Tyz.storeChunk(tyz_data, chunk_offset, chunk_extent);
    Tzx.storeChunk(tzx_data, chunk_offset, chunk_extent);
    Txy.storeChunk(txy_data, chunk_offset, chunk_extent);

    Rho.storeChunk(rho_data, chunk_offset, chunk_extent);
    Ke.storeChunk(ke_data, chunk_offset, chunk_extent);

    series.flush();
  }
};
#endif // VPIC_ENABLE_OPENPMD

/*
   template <typename Policy = BinaryDump>
   struct IODump : private Policy {
   using Policy::dump_particles;
   using Policy::dump_fields;
   using Policy::dump_hydro;
   };
   */

#endif // Dump_Strategy_h
