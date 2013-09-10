#include "common_object.h"
#include <hdf5.h>
#include <hdf5_hl.h>

#define H5FILE_NAME "RUN/field.h5"

void output_h5(){
  char *fname;
  int RANK = 3;
  char NAME[128];
  char NEW_H5FILE_NAME[128];
  char XMF_FILE_NAME[128];
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
  hid_t hdf5_type;
  hsize_t array[ ] = {3}; 

  int i;
  int size = (LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD);
  my_double  *aux;
  
  FILE *fout;

  aux  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
  if(aux == NULL){ fprintf(stderr,"Not enough memory to allocate aux field t\n"); exit(-1);}

  efile_3d[0] = NZ;  efile_3d[1] = NY;  efile_3d[2] = NX;
  efilespace = H5Screate_simple(RANK, efile_3d, NULL);

  edimens_3d[0] = LNZ+TWO_BRD;  edimens_3d[1] = LNY+TWO_BRD;  edimens_3d[2] = LNX+TWO_BRD;
  ememspace = H5Screate_simple(RANK, edimens_3d, NULL);

  plist_id = H5Pcreate(H5P_FILE_ACCESS);

  hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);

  file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  group   = H5Gcreate (file_id, "/euler", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  /*
  file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id>=0) {
    file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, H5P_DEFAULT);
    group   = H5Gopen(file_id, "/euler", H5P_DEFAULT);
  } else {
    file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    group   = H5Gcreate (file_id, "/euler", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  }
  */
  H5Pclose(plist_id);
  
  hdf5_type = H5Tcopy(H5T_NATIVE_DOUBLE);

  property_id  = H5Pcreate(H5P_DATASET_CREATE);       
    
  estart_3d[0] = BRD;  estart_3d[1] = BRD;  estart_3d[2] = BRD;
  estride_3d[0] = 1;  estride_3d[1] = 1;   estride_3d[2] = 1;
  ecount_3d[0] = 1;    ecount_3d[1] = 1;    ecount_3d[2] = 1;
  eblock_3d[0] = LNZ;  eblock_3d[1] = LNY;  eblock_3d[2] = LNX;

  dstart_3d[0] = mez*LNZ;  dstart_3d[1] = mey*LNY;  dstart_3d[2] = mex*LNX;
  dstride_3d[0] = 1;  dstride_3d[1] = 1;  dstride_3d[2] = 1;
  dcount_3d[0] = 1;  dcount_3d[1] = 1;  dcount_3d[2] = 1;
  dblock_3d[0] = LNZ;  dblock_3d[1] = LNY;  dblock_3d[2] = LNX;

  //fprintf(stderr,"mex %d mey %d mez %d\n",mex,mey,mez);

  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d);
  status = H5Sselect_hyperslab(efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
    
  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);

  /* Writing mesh <- the stupid way */
  edataset = H5Dcreate(group, "position_x", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
 for (i=0; i<size; i++) aux[i]=center_V[i].x;
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, aux);
 H5Dclose(edataset);

  edataset = H5Dcreate(group, "position_y", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
 for (i=0; i<size; i++) aux[i]=center_V[i].y;
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, aux);
 H5Dclose(edataset);

 edataset = H5Dcreate(group, "position_z", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
 for (i=0; i<size; i++) aux[i]=center_V[i].z;
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, aux);
 H5Dclose(edataset);

#ifdef LB_FLUID
  edataset = H5Dcreate(group, "velocity_x", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
 for (i=0; i<size; i++) aux[i]=u[i].x;
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, aux);
 H5Dclose(edataset);

  edataset = H5Dcreate(group, "velocity_y", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
 for (i=0; i<size; i++) aux[i]=u[i].y;
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, aux);
 H5Dclose(edataset);

 edataset = H5Dcreate(group, "velocity_z", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
 for (i=0; i<size; i++) aux[i]=u[i].z;
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, aux);
 H5Dclose(edataset);

 edataset = H5Dcreate(group, "density", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, dens);
 H5Dclose(edataset);
#endif

#ifdef LB_TEMPERATURE
 edataset = H5Dcreate(group, "temperature", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, t);
 H5Dclose(edataset);
#endif

#ifdef LB_SCALAR
 edataset = H5Dcreate(group, "scalar", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, s);
 H5Dclose(edataset);
#endif

  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(efilespace);
  H5Sclose(ememspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

  /* we rename the file */

  sprintf(NEW_H5FILE_NAME,"%s/field_%d.h5",OutDir,itime);
  rename (H5FILE_NAME, NEW_H5FILE_NAME);


  /** Part 2: we now write the corresponding XDMF data format **/

  if(ROOT){
  sprintf(XMF_FILE_NAME,"%s/field_%d.xmf",OutDir,itime);
  sprintf(NEW_H5FILE_NAME,"field_%d.h5",itime);
  size=sizeof(my_double);
  fout = fopen(XMF_FILE_NAME,"w");
 
  fprintf(fout,"<?xml version=\"1.0\" ?>\n");
  fprintf(fout,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(fout,"<Xdmf Version=\"2.0\">\n");
  fprintf(fout,"<Domain>\n");
  fprintf(fout,"<Grid Name=\"FiniteVolumeLB\" GridType=\"Uniform\">\n");
  fprintf(fout,"<Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n",NZ,NY,NX);
  fprintf(fout,"<Geometry GeometryType=\"X_Y_Z\">\n");
  fprintf(fout,"<DataItem Name=\"X\" Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/position_x\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"<DataItem Name=\"Y\" Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/position_y\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"<DataItem Name=\"Z\" Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/position_z\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Geometry>\n");

#ifdef LB_FLUID
  fprintf(fout,"<Attribute Name=\"velocity_x\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/velocity_x\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");

  fprintf(fout,"<Attribute Name=\"velocity_y\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/velocity_y\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");

  fprintf(fout,"<Attribute Name=\"velocity_z\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/velocity_z\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");

  fprintf(fout,"<Attribute Name=\"density\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/density\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");
#endif

#ifdef LB_TEMPERATURE
  fprintf(fout,"<Attribute Name=\"temperature\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/temperature\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");
#endif

#ifdef LB_SCALAR
  fprintf(fout,"<Attribute Name=\"scalar\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/scalar\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");
#endif

  fprintf(fout,"<Time Value=\" %e\" />\n",time_now);
  fprintf(fout,"</Grid>\n");
  fprintf(fout,"</Domain>\n");
  fprintf(fout,"</Xdmf>\n");

  fclose(fout);
  }/* end of Xdmf */

}


/******************************************************/

void write_pop_h5(){
  char *fname;
  int RANK = 3;  
  char label[128] , label2[128];
  char NEW_H5FILE_NAME[128];
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
  hid_t hdf5_type;
  hsize_t array[ ] = {3}; 

  int i;
  int size = (LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD);
  my_double  *aux;


  efile_3d[0] = NZ;  efile_3d[1] = NY;  efile_3d[2] = NX;
  efilespace = H5Screate_simple(RANK, efile_3d, NULL);

  edimens_3d[0] = LNZ+TWO_BRD;  edimens_3d[1] = LNY+TWO_BRD;  edimens_3d[2] = LNX+TWO_BRD;
  ememspace = H5Screate_simple(RANK, edimens_3d, NULL);

  plist_id = H5Pcreate(H5P_FILE_ACCESS);

  hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);

  file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  group   = H5Gcreate (file_id, "/population", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  /*
  file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id>=0) {
    file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, H5P_DEFAULT);
    group   = H5Gopen(file_id, "/euler", H5P_DEFAULT);
  } else {
    file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    group   = H5Gcreate (file_id, "/euler", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  }
  */
  H5Pclose(plist_id);
  
  //hdf5_type = H5Tcopy(H5T_NATIVE_DOUBLE);
 hdf5_type = H5Tcreate (H5T_COMPOUND, sizeof(pop));
  for (i = 0; i < NPOP; i++){
    sprintf(label,"p%d",i);
    H5Tinsert(hdf5_type, label, HOFFSET(pop, p[i]), H5T_NATIVE_DOUBLE);
  }


  property_id  = H5Pcreate(H5P_DATASET_CREATE);       
    
  estart_3d[0] = BRD;  estart_3d[1] = BRD;  estart_3d[2] = BRD;
  estride_3d[0] = 1;  estride_3d[1] = 1;   estride_3d[2] = 1;
  ecount_3d[0] = 1;    ecount_3d[1] = 1;    ecount_3d[2] = 1;
  eblock_3d[0] = LNZ;  eblock_3d[1] = LNY;  eblock_3d[2] = LNX;

  dstart_3d[0] = mez*LNZ;  dstart_3d[1] = mey*LNY;  dstart_3d[2] = mex*LNX;
  dstride_3d[0] = 1;  dstride_3d[1] = 1;  dstride_3d[2] = 1;
  dcount_3d[0] = 1;  dcount_3d[1] = 1;  dcount_3d[2] = 1;
  dblock_3d[0] = LNZ;  dblock_3d[1] = LNY;  dblock_3d[2] = LNX;

  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d);
  status = H5Sselect_hyperslab(efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
    
  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);


#ifdef LB_FLUID
  /////edataset = H5Dcreate(file_id, "velocity", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT); 
  edataset = H5Dcreate(group, "velocity", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT); 
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, p);
  H5Dclose(edataset);
#endif

#ifdef LB_TEMPERATURE
  ////  edataset = H5Dcreate(file_id, "temperature", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT); 
  edataset = H5Dcreate(group, "temperature", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT); 
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, g);
  H5Dclose(edataset);
#endif

#ifdef LB_SCALAR
  ////edataset = H5Dcreate(file_id, "scalar", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT); 
  edataset = H5Dcreate(group, "scalar", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT); 
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, h);
  H5Dclose(edataset);
#endif


  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(efilespace);
  H5Sclose(ememspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

  /* we rename the file */

  sprintf(NEW_H5FILE_NAME,"pop.h5");
  rename (H5FILE_NAME, NEW_H5FILE_NAME);

}

/**********************************************/

#define IN_H5FILE_NAME "pop.h5"

void read_pop_h5(){
  char *fname;
  int RANK = 3;  
  char label[128] , label2[128];
  //char NEW_H5FILE_NAME[128];
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
  hid_t hdf5_type;
  hsize_t array[ ] = {3}; 
  int i;

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);

  file_id = H5Fopen(IN_H5FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  group   = H5Gopen(file_id, "/population", H5P_DEFAULT);
  H5Pclose(plist_id);

  hdf5_type = H5Tcreate (H5T_COMPOUND, sizeof(pop));
  for (i = 0; i < NPOP; i++){
    sprintf(label,"p%d",i);
    H5Tinsert(hdf5_type, label, HOFFSET(pop, p[i]), H5T_NATIVE_DOUBLE);
  }
  property_id  = H5Pcreate(H5P_DATASET_CREATE); 

  efile_3d[0] = NZ;  efile_3d[1] = NY;  efile_3d[2] = NX;
  efilespace = H5Screate_simple(RANK, efile_3d, NULL);

  edimens_3d[0] = LNZ+TWO_BRD;  edimens_3d[1] = LNY+TWO_BRD;  edimens_3d[2] = LNX+TWO_BRD;
  ememspace = H5Screate_simple(RANK, edimens_3d, NULL);
    
  estart_3d[0] = BRD;  estart_3d[1] = BRD;  estart_3d[2] = BRD;
  estride_3d[0] = 1;  estride_3d[1] = 1;   estride_3d[2] = 1;
  ecount_3d[0] = 1;    ecount_3d[1] = 1;    ecount_3d[2] = 1;
  eblock_3d[0] = LNZ;  eblock_3d[1] = LNY;  eblock_3d[2] = LNX;

  dstart_3d[0] = mez*LNZ;  dstart_3d[1] = mey*LNY;  dstart_3d[2] = mex*LNX;
  dstride_3d[0] = 1;  dstride_3d[1] = 1;  dstride_3d[2] = 1;
  dcount_3d[0] = 1;  dcount_3d[1] = 1;  dcount_3d[2] = 1;
  dblock_3d[0] = LNZ;  dblock_3d[1] = LNY;  dblock_3d[2] = LNX;

  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d);
  status = H5Sselect_hyperslab(efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
    
  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);

#ifdef LB_FLUID
  edataset = H5Dopen(group,"velocity",H5P_DEFAULT); 
  ret = H5Dread(edataset, hdf5_type, ememspace, efilespace, H5P_DEFAULT, p);
  H5Dclose(edataset);
#endif

#ifdef LB_TEMPERATURE
  edataset = H5Dopen(group, "temperature", H5P_DEFAULT);
  ret = H5Dread(edataset, hdf5_type, ememspace, efilespace, H5P_DEFAULT, g);
  H5Dclose(edataset);
#endif

#ifdef LB_SCALAR
  edataset = H5Dopen(group, "scalar", H5P_DEFAULT);
  ret = H5Dread(edataset, hdf5_type, ememspace, efilespace, H5P_DEFAULT, h);
  H5Dclose(edataset);
#endif


  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(efilespace);
  H5Sclose(ememspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

}




