#include "common_object.h"
#include <hdf5.h>


write_scalar_h5(my_double   *t,  char *fname)
{
  int RANK = 3;
  char H5FILE_NAME[256];
  char DATASETNAME[256];
  hid_t file_id, group, edataset, ememspace, hdf5_status;
  hid_t xfer_plist, ret, property_id, efilespace;
  hsize_t efile_3d[3], edimens_3d[3];
  hsize_t estart_3d[3], ecount_3d[3], estride_3d[3], eblock_3d[3];
  hsize_t dstart_3d[3], dcount_3d[3], dstride_3d[3], dblock_3d[3];
  hid_t plist_id;            /* property list identifier */
  herr_t status;
  H5E_auto_t old_func;
  void *old_client_data;
  hid_t string_type, hdf5_poptype, H5_TYPE_IN_USE;
 

  sprintf(H5FILE_NAME,"%s/%s_%d.h5",OutDir,fname,itime);
  sprintf(DATASETNAME,"%s",fname);

  H5_TYPE_IN_USE = H5T_NATIVE_DOUBLE;

  efile_3d[0] = NX;  efile_3d[1] = NY;  efile_3d[2] = NZ;
  efilespace = H5Screate_simple(3, efile_3d, NULL);
  
  edimens_3d[0] = LNX+TWO_BRD;  edimens_3d[1] = LNY+TWO_BRD;  edimens_3d[2] = LNZ+TWO_BRD;
  ememspace = H5Screate_simple(3, edimens_3d, NULL);

  //sprintf(H5FILE_NAME,"%s/%s.h5",basedir,basename);
  plist_id = H5Pcreate(H5P_FILE_ACCESS);

  hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);

  
  //file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, H5P_DEFAULT);
  //if (file_id>=0) {
    // group   = H5Gopen(file_id, "/EULER",H5P_DEFAULT);
    // } else {
    file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    group   = H5Gcreate (file_id, "/EULER", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    // }

  //H5Eset_auto(H5E_DEFAULT,old_func, old_client_data);
  //write_attributes(group,H5FILE_NAME);
  H5Pclose(plist_id);


  //H5_TYPE_IN_USE = H5T_NATIVE_DOUBLE;
  //hdf5_poptype = H5Tcopy(H5_TYPE_IN_USE);  
  property_id  = H5Pcreate(H5P_DATASET_CREATE);

  edataset = H5Dcreate(group, DATASETNAME, H5T_NATIVE_DOUBLE, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);       
  
  estart_3d[0] = BRD;  estart_3d[1] = BRD;  estart_3d[2] = BRD;
  estride_3d[0] = 1;  estride_3d[1] = 1;  estride_3d[2] = 1;
  eblock_3d[0] = LNX;  eblock_3d[1] = LNY;  eblock_3d[2] = LNZ;
  ecount_3d[0] = 1;  ecount_3d[1] = 1;  ecount_3d[2] = 1;

  dstart_3d[0] = mex*LNX;  dstart_3d[1] = mey*LNY;  dstart_3d[2] = mez*LNZ;
  dstride_3d[0] = 1;  dstride_3d[1] = 1;  dstride_3d[2] = 1;
  dblock_3d[0] = LNX;  dblock_3d[1] = LNY;  dblock_3d[2] = LNZ;
  dcount_3d[0] = 1;  dcount_3d[1] = 1;  dcount_3d[2] = 1;

  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d);
  status = H5Sselect_hyperslab(efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
  
  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);
  ret = H5Dwrite(edataset, H5T_NATIVE_DOUBLE, ememspace, efilespace, xfer_plist, t);

  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(efilespace);
  H5Sclose(ememspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  /* Close the data set */
  H5Dclose(edataset);
  H5Gclose(group);
  H5Fclose(file_id);  
}
  


write_vector_h5(vector   *t, char *fname){
  int RANK = 3;
  char H5FILE_NAME[256];
  char DATASETNAME[256];
  hid_t file_id, group, edataset, ememspace, hdf5_status;
  hid_t xfer_plist, ret, property_id, efilespace;
  hsize_t efile_3d[3], edimens_3d[3];
  hsize_t estart_3d[3], ecount_3d[3], estride_3d[3], eblock_3d[3];
  hsize_t dstart_3d[3], dcount_3d[3], dstride_3d[3], dblock_3d[3];
  hid_t plist_id;            /* property list identifier */
  herr_t status;
  H5E_auto_t old_func;
  void *old_client_data;
  hid_t hdf5_pop_type;
 
  sprintf(H5FILE_NAME,"%s/%s_%d.h5",OutDir,fname,itime);
  sprintf(DATASETNAME,"%s",fname);

  //H5_TYPE_IN_USE = H5T_NATIVE_DOUBLE;

  efile_3d[0] = NX;  efile_3d[1] = NY;  efile_3d[2] = NZ;
  efilespace = H5Screate_simple(3, efile_3d, NULL);
  
  edimens_3d[0] = LNX+TWO_BRD;  edimens_3d[1] = LNY+TWO_BRD;  edimens_3d[2] = LNZ+TWO_BRD;
  ememspace = H5Screate_simple(3, edimens_3d, NULL);

  //sprintf(H5FILE_NAME,"%s/%s.h5",basedir,basename);
  plist_id = H5Pcreate(H5P_FILE_ACCESS);

  hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);

  
  //file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, H5P_DEFAULT);
  //if (file_id>=0) {
    // group   = H5Gopen(file_id, "/EULER",H5P_DEFAULT);
    // } else {
    file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    group   = H5Gcreate (file_id, "/EULER", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    // }

  //H5Eset_auto(H5E_DEFAULT,old_func, old_client_data);
  //write_attributes(group,H5FILE_NAME);
  H5Pclose(plist_id);


  //H5_TYPE_IN_USE = H5T_NATIVE_DOUBLE;
  //hdf5_poptype = H5Tcopy(H5_TYPE_IN_USE);  
  hdf5_pop_type = H5Tcreate (H5T_COMPOUND, sizeof(vector));
  H5Tinsert(hdf5_pop_type, "x", HOFFSET(vector, x), H5T_NATIVE_DOUBLE);
  H5Tinsert(hdf5_pop_type, "y", HOFFSET(vector, y), H5T_NATIVE_DOUBLE);
  H5Tinsert(hdf5_pop_type, "z", HOFFSET(vector, z), H5T_NATIVE_DOUBLE);

  property_id  = H5Pcreate(H5P_DATASET_CREATE);

  edataset = H5Dcreate(group, DATASETNAME, hdf5_pop_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);       
  
  estart_3d[0] = BRD;  estart_3d[1] = BRD;  estart_3d[2] = BRD;
  estride_3d[0] = 1;  estride_3d[1] = 1;  estride_3d[2] = 1;
  eblock_3d[0] = LNX;  eblock_3d[1] = LNY;  eblock_3d[2] = LNZ;
  ecount_3d[0] = 1;  ecount_3d[1] = 1;  ecount_3d[2] = 1;

  dstart_3d[0] = mex*LNX;  dstart_3d[1] = mey*LNY;  dstart_3d[2] = mez*LNZ;
  dstride_3d[0] = 1;  dstride_3d[1] = 1;  dstride_3d[2] = 1;
  dblock_3d[0] = LNX;  dblock_3d[1] = LNY;  dblock_3d[2] = LNZ;
  dcount_3d[0] = 1;  dcount_3d[1] = 1;  dcount_3d[2] = 1;

  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d);
  status = H5Sselect_hyperslab(efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
  
  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);
  ret = H5Dwrite(edataset, hdf5_pop_type, ememspace, efilespace, xfer_plist, t);

  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(efilespace);
  H5Sclose(ememspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  /* Close the data set */
  H5Dclose(edataset);
  H5Gclose(group);
  H5Fclose(file_id);  
}




