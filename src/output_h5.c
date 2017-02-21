#include "common_object.h"
#include <hdf5.h>
#include <hdf5_hl.h>

//#define H5FILE_NAME "RUN/field.h5"

#define ARANK  1   /* Rank and dimension sizes of the first dataset attribute */

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

  // hsize_t     family_size = 2; // family
  //#define FAMIGLIA "%d"
  //sprintf(NEW_H5FILE_NAME,"%s/field_%d_%s.h5",OutDir,itime,FAMIGLIA); //family
  sprintf(NEW_H5FILE_NAME,"%s/field_%d.h5",OutDir,itime);
  if(ROOT) fprintf(stderr,"Writing file %s\n",NEW_H5FILE_NAME);


  aux  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
  if(aux == NULL){ fprintf(stderr,"Not enough memory to allocate aux field t\n"); exit(-1);}

  efile_3d[0] = NZ;  efile_3d[1] = NY;  efile_3d[2] = NX;
  efilespace = H5Screate_simple(RANK, efile_3d, NULL);

  edimens_3d[0] = LNZ+TWO_BRD;  edimens_3d[1] = LNY+TWO_BRD;  edimens_3d[2] = LNX+TWO_BRD;
  ememspace = H5Screate_simple(RANK, edimens_3d, NULL);

  plist_id = H5Pcreate(H5P_FILE_ACCESS);

  //status = H5Pset_fapl_family(plist_id, family_size, H5P_DEFAULT);//family
  //if(ROOT) printf ("H5Pset_fapl_family returns: %i\n", status); // family

  hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);

  file_id = H5Fcreate(NEW_H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
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

#ifdef OUTPUT_H5_GRID
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
#endif

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
 #ifdef LB_TEMPERATURE_MELTING
  edataset = H5Dcreate(group, "liquid_fraction", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, liquid_frac);
  H5Dclose(edataset);
 #endif
#endif

#ifdef LB_SCALAR
 edataset = H5Dcreate(group, "scalar", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, s);
 H5Dclose(edataset);
#endif

#ifdef LB_FLUID_FORCING_LANDSCAPE
 edataset = H5Dcreate(group, "landscape", hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
  ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, landscape);
 H5Dclose(edataset);
#endif



  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(efilespace);
  H5Sclose(ememspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

  /* free scalar auxiliary field */
  free(aux);

  /* we rename the file */

  //sprintf(NEW_H5FILE_NAME,"%s/field_%d.h5",OutDir,itime);
  //if(ROOT) rename(H5FILE_NAME, NEW_H5FILE_NAME);


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
#ifdef OUTPUT_H5_GRID
  /* we write the grid coordinates - useful for the finite volume approach */ 
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
#else
  fprintf(fout,"<Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",NZ,NY,NX);
  fprintf(fout,"<Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
  fprintf(fout,"<DataItem Name=\"Z_Y_X\" Dimensions=\"3\" DataType=\"Float\" Format=\"XML\">\n");
  fprintf(fout,"0.5 0.5 0.5\n");
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"<DataItem Name=\"DZ_DY_DX\" Dimensions=\"3\" DataType=\"Float\" Format=\"XML\">\n");
  fprintf(fout,"1.0 1.0 1.0\n");
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Geometry>\n");
#endif

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

  /* to read the velocity as a vector */
  fprintf(fout,"<Attribute Name=\"velocity\" AttributeType=\"Vector\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem ItemType=\"Function\" Dimensions=\"%d %d %d 3\" \n   Function=\"JOIN($0 , $1, $2)\">\n",NZ,NY,NX);
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/velocity_x\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/velocity_y\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d 1\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/velocity_z\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");
#endif

#ifdef LB_TEMPERATURE
  fprintf(fout,"<Attribute Name=\"temperature\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/temperature\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");
 #ifdef LB_TEMPERATURE_MELTING
  fprintf(fout,"<Attribute Name=\"liquid_fraction\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/liquid_fraction\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");
 #endif
#endif

#ifdef LB_SCALAR
  fprintf(fout,"<Attribute Name=\"scalar\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/scalar\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");
#endif

#ifdef LB_FLUID_FORCING_LANDSCAPE
  fprintf(fout,"<Attribute Name=\"landscape\" AttributeType=\"Scalar\" Center=\"Node\">\n");
  fprintf(fout,"<DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",NZ,NY,NX,size);
  fprintf(fout,"%s:/euler/landscape\n",NEW_H5FILE_NAME);
  fprintf(fout,"</DataItem>\n");
  fprintf(fout,"</Attribute>\n");
#endif

  fprintf(fout,"<Time Value=\" %e\" />\n",time_now);
  fprintf(fout,"</Grid>\n");
  fprintf(fout,"</Domain>\n");
  fprintf(fout,"</Xdmf>\n");

  fclose(fout);
  }/* end of Xdmf */

  MPI_Barrier(MPI_COMM_WORLD);

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
  //my_double  *aux;

  /* definitions for attributes */
#ifdef LB_FLUID_FORCING_HIT
   hid_t   attr1, attr2, attr3; /* Attribute identifiers */
   hid_t   attr;
   hid_t   aid1, aid2, aid3;    /* Attribute dataspace identifiers */ 
   hid_t   atype;               /* Attribute type */

   //hsize_t fdim[] = {SIZE};
   hsize_t adim[] = {nk};  /* Dimensions of the first attribute  */
#endif

  sprintf(NEW_H5FILE_NAME,"pop.h5");
  if(ROOT) fprintf(stderr,"Writing file %s\n",NEW_H5FILE_NAME);

  efile_3d[0] = NZ;  efile_3d[1] = NY;  efile_3d[2] = NX;
  efilespace = H5Screate_simple(RANK, efile_3d, NULL);

  edimens_3d[0] = LNZ+TWO_BRD;  edimens_3d[1] = LNY+TWO_BRD;  edimens_3d[2] = LNX+TWO_BRD;
  ememspace = H5Screate_simple(RANK, edimens_3d, NULL);

  plist_id = H5Pcreate(H5P_FILE_ACCESS);

  hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);

  file_id = H5Fcreate(NEW_H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
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

/* ATTRIBUTES for FORCING */
#if defined(LB_FLUID_FORCING_HIT)||defined(LB_TEMPERATURE_FORCING_HIT)||defined(LB_SCALAR_FORCING_HIT) 

 hdf5_type = H5Tcreate (H5T_COMPOUND, sizeof(vector));
 H5Tinsert(hdf5_type, "x", HOFFSET(vector, x), H5T_NATIVE_DOUBLE);
 H5Tinsert(hdf5_type, "y", HOFFSET(vector, y), H5T_NATIVE_DOUBLE);
 H5Tinsert(hdf5_type, "z", HOFFSET(vector, z), H5T_NATIVE_DOUBLE);

 /* Create a group */
 group   = H5Gcreate (file_id, "/forcing", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

 /* Create scalar attribute : number of modes */
 aid1  = H5Screate(H5S_SCALAR);
 attr1 = H5Acreate(group, "number of modes (nk)", H5T_NATIVE_INT, aid1,H5P_DEFAULT,H5P_DEFAULT);
 /* Write scalar attribute */
 ret = H5Awrite(attr1, H5T_NATIVE_INT, &nk); 
 H5Aclose(attr1);
 H5Sclose(aid1); 
 

 #ifdef LB_FLUID_FORCING_HIT
   /* Create dataspace for the phase attirbutes */
   aid2 = H5Screate(H5S_SIMPLE);
   ret  = H5Sset_extent_simple(aid2, ARANK, adim, NULL);
   /* Create array attribute */
   attr2 = H5Acreate(group, "phi_u", hdf5_type, aid2,H5P_DEFAULT,H5P_DEFAULT);
   /* Write array attribute */
   ret = H5Awrite(attr2, hdf5_type, phi_u);
   H5Aclose(attr2);
   H5Sclose(aid2);
 #endif

 #ifdef LB_TEMPERATURE_FORCING_HIT
   /* Create dataspace for the phase attirbutes */
   aid2 = H5Screate(H5S_SIMPLE);
   ret  = H5Sset_extent_simple(aid2, ARANK, adim, NULL);
   /* Create array attribute */
   attr2 = H5Acreate(group, "phi_t", hdf5_type, aid2,H5P_DEFAULT,H5P_DEFAULT);
   /* Write array attribute */
   ret = H5Awrite(attr2, hdf5_type, phi_t);
   H5Aclose(attr2);
   H5Sclose(aid2);
 #endif

 #ifdef LB_SCALAR_FORCING_HIT
   /* Create dataspace for the phase attirbutes */
   aid2 = H5Screate(H5S_SIMPLE);
   ret  = H5Sset_extent_simple(aid2, ARANK, adim, NULL);
   /* Create array attribute */
   attr2 = H5Acreate(group, "phi_s", hdf5_type, aid2,H5P_DEFAULT,H5P_DEFAULT);
   /* Write array attribute */
   ret = H5Awrite(attr2, hdf5_type, phi_s);
   H5Aclose(attr2);
   H5Sclose(aid2);
 #endif
  
   /* close group */
   H5Gclose(group);
#endif


  H5Fclose(file_id);  

  /* we rename the file */

  //  sprintf(NEW_H5FILE_NAME,"pop.h5");
  //  rename (H5FILE_NAME, NEW_H5FILE_NAME);

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

#ifdef LB_FLUID_FORCING_HIT
   hid_t   attr1, attr2, attr3; /* Attribute identifiers */
   hid_t   attr;
   hid_t   aid1, aid2, aid3;    /* Attribute dataspace identifiers */ 
   hid_t   atype , aspace;               /* Attribute type */
   

   //hsize_t fdim[] = {SIZE};
   hsize_t adim[] = {nk};  /* Dimensions of the first attribute  */
#endif

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
  if(resume_u){
  edataset = H5Dopen(group,"velocity",H5P_DEFAULT); 
  ret = H5Dread(edataset, hdf5_type, ememspace, efilespace, H5P_DEFAULT, p);
  H5Dclose(edataset);
  }
#endif

#ifdef LB_TEMPERATURE
  if(resume_t){
  edataset = H5Dopen(group, "temperature", H5P_DEFAULT);
  ret = H5Dread(edataset, hdf5_type, ememspace, efilespace, H5P_DEFAULT, g);
  H5Dclose(edataset);
  }
#endif

#ifdef LB_SCALAR
  if(resume_s){
  edataset = H5Dopen(group, "scalar", H5P_DEFAULT);
  ret = H5Dread(edataset, hdf5_type, ememspace, efilespace, H5P_DEFAULT, h);
  H5Dclose(edataset);
  }
#endif

  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(efilespace);
  H5Sclose(ememspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);


/* ATTRIBUTES for FORCING */
#if defined(LB_FLUID_FORCING_HIT)||defined(LB_TEMPERATURE_FORCING_HIT)||defined(LB_SCALAR_FORCING_HIT) 
   hdf5_type = H5Tcreate (H5T_COMPOUND, sizeof(vector));
   H5Tinsert(hdf5_type, "x", HOFFSET(vector, x), H5T_NATIVE_DOUBLE);
   H5Tinsert(hdf5_type, "y", HOFFSET(vector, y), H5T_NATIVE_DOUBLE);
   H5Tinsert(hdf5_type, "z", HOFFSET(vector, z), H5T_NATIVE_DOUBLE);

/* Check existence of the group "/forcing"  */
   status = H5Eset_auto(NULL, NULL,NULL);
   status = H5Gget_objinfo (file_id, "/forcing", 0, NULL);
   if (status == 0){ 
    if(ROOT) printf ("The group /forcing exists.\n");

 /* a group */
 group   = H5Gopen(file_id, "/forcing", H5P_DEFAULT);

 property_id  = H5Pcreate(H5P_GROUP_ACCESS);
 //xfer_plist = H5Pcreate(H5P_DATASET_XFER);
 //ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE); 

 /* Create attirbute identifier */
 attr1 = H5Aopen(group, "number of modes (nk)",H5P_DEFAULT);
 /* Read scalar attribute */
 ret = H5Aread(attr1, H5T_NATIVE_INT, &nk); 
 H5Aclose(attr1);
 if(ROOT)fprintf(stderr,"The value of the attribute \"number of modes (nk)\" is %d \n", nk);

 #ifdef LB_FLUID_FORCING_HIT
  if(resume_u){ 
  /* Create attirbute identifier */
   attr2 = H5Aopen(group, "phi_u", H5P_DEFAULT);
   //aspace = H5Aget_space(attr2);
   //atype  = H5Aget_type(attr2);
   phi_u  = (vector*) malloc(sizeof(vector)*nk);
   /* Read array attribute */
   ret = H5Aread(attr2,hdf5_type, phi_u);
   //H5Tclose(atype);
   H5Aclose(attr2);
   if(ROOT) for (i = 0; i < nk; i++) fprintf(stderr,"The value of the attribute \"phi_u\" is %e %e %e \n", phi_u[i].x,phi_u[i].y ,phi_u[i].z);
  }
 #endif

 #ifdef LB_TEMPERATURE_FORCING_HIT
  if(resume_t){ 
   /* Create attirbute identifier */
   attr2 = H5Aopen(group, "phi_t", H5P_DEFAULT);
   //aspace = H5Aget_space(attr2);
   //atype  = H5Aget_type(attr2);
   phi_t  = (vector*) malloc(sizeof(vector)*nk);
   /* Read array attribute */
   ret = H5Aread(attr2,hdf5_type, phi_t);
   //H5Tclose(atype);
   H5Aclose(attr2);
   if(ROOT) for (i = 0; i < nk; i++) fprintf(stderr,"The value of the attribute \"phi_t\" is %e %e %e \n", phi_t[i].x,phi_t[i].y ,phi_t[i].z);
  }
 #endif

 #ifdef LB_SCALAR_FORCING_HIT
  if(resume_s){ 
    /* Create attirbute identifier */
   attr2 = H5Aopen(group, "phi_s", H5P_DEFAULT);
   //aspace = H5Aget_space(attr2);
   //atype  = H5Aget_type(attr2);
   phi_s  = (vector*) malloc(sizeof(vector)*nk);
   /* Read array attribute */
   ret = H5Aread(attr2,hdf5_type, phi_s);
   //H5Tclose(atype);
   H5Aclose(attr2);
   if(ROOT) for (i = 0; i < nk; i++) fprintf(stderr,"The value of the attribute \"phi_s\" is %e %e %e \n", phi_s[i].x,phi_s[i].y ,phi_s[i].z);
  }
 #endif
  
   /* close group */
   H5Pclose(property_id);
   H5Gclose(group);

 }else{
  if(ROOT)  fprintf(stderr,"HDF5: Group \"/forcing\" not found! \n");
 }
 #endif

  /* close the file */
  H5Fclose(file_id);  

}

/*******************************************************************************/



void write_scalar_h5(my_double *field, char which_field[128]){
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

  sprintf(NEW_H5FILE_NAME,"%s.h5", which_field);
  if(ROOT) fprintf(stderr,"Writing file %s\n",NEW_H5FILE_NAME);

  aux  = (my_double*) malloc(sizeof(my_double)*(LNX+TWO_BRD)*(LNY+TWO_BRD)*(LNZ+TWO_BRD)); 
  if(aux == NULL){ fprintf(stderr,"Not enough memory to allocate aux field t\n"); exit(-1);}

  efile_3d[0] = NZ;  efile_3d[1] = NY;  efile_3d[2] = NX;
  efilespace = H5Screate_simple(RANK, efile_3d, NULL);

  edimens_3d[0] = LNZ+TWO_BRD;  edimens_3d[1] = LNY+TWO_BRD;  edimens_3d[2] = LNX+TWO_BRD;
  ememspace = H5Screate_simple(RANK, edimens_3d, NULL);

  plist_id = H5Pcreate(H5P_FILE_ACCESS);

  hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);

  file_id = H5Fcreate(NEW_H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  group   = H5Gcreate (file_id, "/euler", H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
 
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

  status = H5Sselect_hyperslab( ememspace, H5S_SELECT_SET, estart_3d, estride_3d, ecount_3d, eblock_3d);
  status = H5Sselect_hyperslab(efilespace, H5S_SELECT_SET, dstart_3d, dstride_3d, dcount_3d, dblock_3d);
    
  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);

 /* Here we write just the field */
 edataset = H5Dcreate(group, which_field, hdf5_type, efilespace,H5P_DEFAULT, property_id,H5P_DEFAULT);
 ret = H5Dwrite(edataset, hdf5_type, ememspace, efilespace, xfer_plist, field);
 H5Dclose(edataset);


  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(efilespace);
  H5Sclose(ememspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

  /* free scalar auxiliary field */
  free(aux);

  /* we rename the file */

  //  sprintf(NEW_H5FILE_NAME,"%s.h5", which_field);
  //  if(ROOT) rename(H5FILE_NAME, NEW_H5FILE_NAME);

  MPI_Barrier(MPI_COMM_WORLD);
}


/******************************************************/

void read_scalar_h5(my_double *field, char which_field[128]){
  char *fname;
  int RANK = 3;  
  char label[128] , label2[128];
  char SCALAR_H5FILE_NAME[128];
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
  FILE *fin;

  sprintf(SCALAR_H5FILE_NAME,"%s.h5", which_field);

  if(ROOT){
    fin = fopen(SCALAR_H5FILE_NAME,"r");
    if(fin == NULL){
	 fprintf(stderr,"Error message -> %s file is missing! Exit.\n",SCALAR_H5FILE_NAME);
	 fclose(fin);
	 exit(0);
    }
    fclose(fin);
    fprintf(stderr,"INITIAL CONDITION : reading file  %s\n",SCALAR_H5FILE_NAME);
  }

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  hdf5_status  = H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD,  MPI_INFO_NULL);

  file_id = H5Fopen(SCALAR_H5FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  group   = H5Gopen(file_id, "/euler", H5P_DEFAULT);
  H5Pclose(plist_id);

  hdf5_type = H5Tcopy(H5T_NATIVE_DOUBLE);

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

  edataset = H5Dopen(group, which_field, H5P_DEFAULT);
  ret = H5Dread(edataset, hdf5_type, ememspace, efilespace, H5P_DEFAULT, field);
  H5Dclose(edataset);


  MPI_Barrier(MPI_COMM_WORLD);
      
  H5Sclose(efilespace);
  H5Sclose(ememspace);
  H5Pclose(xfer_plist);
  H5Pclose(property_id);
  H5Gclose(group);
  H5Fclose(file_id);  

}

/******************************/
