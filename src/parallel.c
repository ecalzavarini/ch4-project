#include "common_object.h"


void initialization_MPI(int *argc, char ***argv){
	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if (ROOT)
		fprintf(stderr, "Hello world! I am ROOT: processor %d of %d\n", me, nprocs);


	/* commit types */
MPI_Type_contiguous(NPOP, MPI_DOUBLE, &MPI_pop_type);
MPI_Type_commit(&MPI_pop_type);

MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_vector_type);
MPI_Type_commit(&MPI_vector_type);
}


/**********************************************************************************/

int 
compare(const void *a, const void *b)
{
	return (*(int *) a - *(int *) b);
}



void 
processor_splitting()
{

	int             tnprocs, nleft, ncubic;
	int             sort3d[3], sort2d[2];
	int             i, j, k, ami;

	/* grid processor splitting */
	LNX = NX;
	LNY = NY;
	LNZ = NZ;
	tnprocs = nprocs;
	nxprocs = nyprocs = nzprocs = 1;

	if (ROOT) {
		while (nxprocs * nyprocs * nzprocs < nprocs) {
			sort3d[0] = (LNX % 2 == 0) ? LNX : 0;
			sort3d[1] = (LNY % 2 == 0) ? LNY : 0;
			sort3d[2] = (LNZ % 2 == 0) ? LNZ : 0;
			qsort(sort3d, 3, sizeof(int), compare);
			if (sort3d[2] == LNX) {
				LNX /= 2;
				nxprocs *= 2;
			} else {
				if (sort3d[2] == LNY) {
					LNY /= 2;
					nyprocs *= 2;
				} else {
					if (sort3d[2] == LNZ) {
						LNZ /= 2;
						nzprocs *= 2;
					}
				}
			}
			if (LNX % 2 != 0 && LNY % 2 != 0 && LNZ % 2 != 0)
				break;
		}
		if (nxprocs * nyprocs * nzprocs == nprocs) {
			if (!me)
				fprintf(stderr, "good splitting!\n");
		} else {
			if (!me)
				fprintf(stderr, "bad splitting: %d over %d\n", nxprocs * nyprocs * nzprocs, nprocs);
		}
		if (!me)
			fprintf(stderr, "LNX %d LNY %d LNZ %d\n", LNX, LNY, LNZ);
		if (!me)
			fprintf(stderr, "nxprocs %d nyprocs %d nzprocs %d total %d\n", nxprocs, nyprocs, nzprocs, nxprocs * nyprocs * nzprocs);
	}			/* if ROOT */
	/* Now broadcast */
	MPI_Bcast(&LNX, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&LNY, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&LNZ, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nxprocs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nyprocs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nzprocs, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef DEBUG
	for (i = 0; i < nprocs; i++) {
		if (i == me)
			fprintf(stderr, "me %d , System Size: LNX %d LNY %d LNZ %d\n   nxprocs %d nyprocs %d nzprocs %d\n", me, LNX, LNY, LNZ, nxprocs, nyprocs, nzprocs);
	}
#endif

	/* here we find the coordinate of the processors */
	for (k = 0; k < nzprocs; k++)
		for (j = 0; j < nyprocs; j++)
			for (i = 0; i < nxprocs; i++) {
				ami = k * (nyprocs * nxprocs) + j * nxprocs + i;
				if (ami == me) {
					mex = i;
					mey = j;
					mez = k;
				}
			}

#ifdef DEBUG
	fprintf(stderr, "me %d , mex %d  mey %d mez %d\n", me, mex, mey, mez);
#endif

	/* processor rulers for vertices */
	LNX_START = LNX * mex;
	LNX_END = LNX * (mex + 1);
	LNY_START = LNY * mey;
	LNY_END = LNY * (mey + 1);
	LNZ_START = LNZ * mez;
	LNZ_END = LNZ * (mez + 1);

	/* processor rulers for the grid */
	NXG=NX+1;
	NYG=NY+1;
	NZG=NZ+1;
	LNXG=LNX+1;
	LNYG=LNY+1;
	LNZG=LNZ+1;
	LNXG_START = LNXG * mex;
	LNXG_END = LNXG * (mex + 1);
	LNYG_START = LNYG * mey;
	LNYG_END = LNYG * (mey + 1);
	LNZG_START = LNZG * mez;
	LNZG_END = LNZG * (mez + 1);

	/* every procs finds it neighbors */
	me_xp = (mex+1+nxprocs) % nxprocs;  
	me_xm = (mex-1+nxprocs) % nxprocs;
	me_yp = (mey+1+nyprocs) % nyprocs;  
	me_ym = (mey-1+nyprocs) % nyprocs;
	me_zp = (mez+1+nzprocs) % nzprocs;  
	me_zm = (mez-1+nzprocs) % nzprocs;

}



