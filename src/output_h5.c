#include "common_object.h"
#include <hdf5.h>


#define H5FILE_NAME     "temp.h5"
#define DATASETNAME 	"temperature"
#define RANK   3       /* dataset rank */


void write_h5(my_double   *t){
 

    /*
     * HDF5 APIs definitions
     */ 	
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[] = {NX, NY,NZ};                 /* dataset dimensions */
    hid_t	plist_id;                 /* property list identifier */
    int         i;
    herr_t      status;

    /* 
     * Set up file access property list with parallel I/O access
     */
     plist_id = H5Pcreate(H5P_FILE_ACCESS);
     H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    /*
     * Create a new file collectively and release property list identifier.
     */
    file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
              H5Pclose(plist_id);
   
    /*
     * Create the dataspace for the dataset.
     */
    filespace = H5Screate_simple(RANK, dimsf, NULL); 

    /*
     * Create the dataset with default properties and close filespace.
     */
    dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_DOUBLE, filespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    /*
     * Create property list for collective dataset write.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
    /*
     * To write dataset independently use
     *
     * H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); 
     */
    
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id, t);
    /*
    free(data);
    */
    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Pclose(plist_id);
    H5Fclose(file_id);
 
}
  


  







