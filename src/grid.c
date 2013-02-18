#include "common_object.h"

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

	/* processor rulers */
	LNX_START = LNX * mex;
	LNX_END = LNX * (mex + 1);
	LNY_START = LNY * mey;
	LNY_END = LNY * (mey + 1);
	LNZ_START = LNZ * mez;
	LNZ_END = LNZ * (mez + 1);

}



void read_mesh()
{

	char            fnamein[256], fnameout[256];
	char            name[256] = "NULL";
	FILE           *fin, *fout;
	int             i, j, k;

	sprintf(fnamein, "mesh.in");
	fin = fopen(fnamein, "r");
	if (fin != NULL) {
		fprintf(stderr, "Mesh file %s has been found!\n", fnamein);
	} else {

		fprintf(stderr, "Warning message -> %s file is missing!\n Starting from grid generated on the fly\n ", fnamein);

		/* moving on the bulk only */
		for (k = 1; k < (LNZ + 1)+BZ-1; k++)
			for (j = 1; j < (LNY+1)+BY-1; j++)
				for (i = 1; i < (LNX + 1)+BX-1; i++) {
					mesh[IDXG(i, j, k)].x = (my_double) (i + LNX_START- BX/2);
					mesh[IDXG(i, j, k)].y = (my_double) (j + LNY_START- BY/2);
					mesh[IDXG(i, j, k)].z = (my_double) (k + LNZ_START- BZ/2);
					/*
					 * flag: 1 is bulk , 0 is wall , -1
					 * is dormient
					 */
					mesh[IDXG(i, j, k)].flag = 1;
				}	/* for ijk */

	}			/* endif */


#ifdef DEBUG
	/* Each processor prints its mesh */
	sprintf(fnamein, "mesh.%d.out", me);
	fout = fopen(fnamein, "w");
	for (k = 0; k < LNZ + 1 + BX; k++)
		for (j = 0; j < LNY+ 1 + BY; j++)
			for (i = 0; i < LNX + 1 + BZ; i++)
				fprintf(fout, "%d %d %d %e %e %e\n", i, j, k, mesh[IDXG(i, j, k)].x, mesh[IDXG(i, j, k)].y, mesh[IDXG(i, j, k)].z);
	fclose(fout);
#endif

}



void 
compute_volumes()
{
	int             i, j, k, pp;

	/* Allocating array for storing information of control volume */
	my_double       P0[3], P1[3], P2[3], P3[3], P4[3], P5[3], P6[3],
	                P7[3], P8[3];
	my_double       D17[3], D35[3], D03[3], D12[3], D05[3], D14[3],
	                D06[3], D24[3], D27[3], D36[3], D47[3], D56[3];
	int             n;
	my_double       S1357, S0145, S0246, S0123, S4567, S2367;
	my_double       N1357[3], N0145[3], N0246[3], N0123[3], N4567[3],
	                N2367[3];
	my_double       M0167[3], M0257[3], M0347[3];
	my_double       V, V1, V2, V3;



	for (i = 0; i < LNX + BX ; i++) {
		for (j = 0; j < LNY + BY; j++) {
			for (k = 0; k < LNZ + BZ; k++) {

/* points definition */
/*    
          z   
          ^ 
          |                         4 -----6
          |                        /|    / |
          | ------>y              / |   /  | 
         /                       5-----7   |
        /                        |  0--|---2 
       /                         | /   |  /
      /                          |/    |/
     x                           1-----3

 */

				P0[0] = mesh[IDXG(i, j, k)].x;
				P0[1] = mesh[IDXG(i, j, k)].y;
				P0[2] = mesh[IDXG(i, j, k)].z;
				P1[0] = mesh[IDXG(i + 1, j, k)].x;
				P1[1] = mesh[IDXG(i + 1, j, k)].y;
				P1[2] = mesh[IDXG(i + 1, j, k)].z;
				P2[0] = mesh[IDXG(i, j + 1, k)].x;
				P2[1] = mesh[IDXG(i, j + 1, k)].y;
				P2[2] = mesh[IDXG(i, j + 1, k)].z;
				P3[0] = mesh[IDXG(i + 1, j + 1, k)].x;
				P3[1] = mesh[IDXG(i + 1, j + 1, k)].y;
				P3[2] = mesh[IDXG(i + 1, j + 1, k)].z;
				P4[0] = mesh[IDXG(i, j, k + 1)].x;
				P4[1] = mesh[IDXG(i, j, k + 1)].y;
				P4[2] = mesh[IDXG(i, j, k + 1)].z;
				P5[0] = mesh[IDXG(i + 1, j, k + 1)].x;
				P5[1] = mesh[IDXG(i + 1, j, k + 1)].y;
				P5[2] = mesh[IDXG(i + 1, j, k + 1)].z;
				P6[0] = mesh[IDXG(i, j + 1, k + 1)].x;
				P6[1] = mesh[IDXG(i, j + 1, k + 1)].y;
				P6[2] = mesh[IDXG(i, j + 1, k + 1)].z;
				P7[0] = mesh[IDXG(i + 1, j + 1, k + 1)].x;
				P7[1] = mesh[IDXG(i + 1, j + 1, k + 1)].y;
				P7[2] = mesh[IDXG(i + 1, j + 1, k + 1)].z;
				P8[0] = (P0[0] + P1[0] + P2[0] + P3[0] + P4[0] + P5[0] + P6[0] + P7[0]) / 8;
				P8[1] = (P0[1] + P1[1] + P2[1] + P3[1] + P4[1] + P5[1] + P6[1] + P7[1]) / 8;
				P8[2] = (P0[2] + P1[2] + P2[2] + P3[2] + P4[2] + P5[2] + P6[2] + P7[2]) / 8;
				for (n = 0; n < 3; n++) {
					fprintf(stdout, "%e %e %e %e %e %e %e %e %e \n", P0[n], P1[n], P2[n], P3[n], P4[n], P5[n], P6[n], P7[n], P8[n]);
				}
				/*
				 * diagonals definition of each surface of a
				 * polyhedron treating as a vector
				 */
				for (n = 0; n < 3; n++) {
					D17[n] = P7[n] - P1[n];
					D35[n] = P5[n] - P3[n];
					D03[n] = P3[n] - P0[n];
					D12[n] = P2[n] - P1[n];
					D05[n] = P5[n] - P0[n];
					D14[n] = P4[n] - P1[n];
					D06[n] = P6[n] - P0[n];
					D24[n] = P4[n] - P2[n];
					D27[n] = P7[n] - P2[n];
					D36[n] = P6[n] - P3[n];
					D47[n] = P7[n] - P4[n];
					D56[n] = P6[n] - P5[n];
					fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e \n", D17[n], D35[n], D03[n], D12[n], D05[n], D14[n], D06[n], D24[n], D27[n], D36[n], D47[n], D56[n]);
				}


				/*
				 * SURFACE AREA definition  Vector formulas
				 * for a quadrilateral: The area of a
				 * quadrilateral ABCD can be calculated using
				 * vectors. Let vectors AC and BD form the
				 * diagonals from A to C and from B to D. The
				 * area of the quadrilateral is then S = 0.5 *
				 * mod(AC X BD) and the magnitude of AC X BD
				 * is calculated by finding the determinant
				 * of matrix formed
				 */
				S1357 = 0.5 * fabs(sqrt((pow((D17[1] * D35[2] - D17[2] * D35[1]), 2)) + (pow((D17[2] * D35[0] - D17[0] * D35[2]), 2)) + (pow((D17[0] * D35[1] - D17[1] * D35[0]), 2))));
				S0145 = 0.5 * fabs(sqrt((pow((D05[1] * D14[2] - D05[2] * D14[1]), 2)) + (pow((D05[2] * D14[0] - D05[0] * D14[2]), 2)) + (pow((D05[0] * D14[1] - D05[1] * D14[0]), 2))));
				S0246 = 0.5 * fabs(sqrt((pow((D06[1] * D24[2] - D06[2] * D24[1]), 2)) + (pow((D06[2] * D24[0] - D06[0] * D24[2]), 2)) + (pow((D06[0] * D24[1] - D06[1] * D24[0]), 2))));
				S2367 = 0.5 * fabs(sqrt((pow((D27[1] * D36[2] - D27[2] * D36[1]), 2)) + (pow((D27[2] * D36[0] - D27[0] * D36[2]), 2)) + (pow((D27[0] * D36[1] - D27[1] * D36[0]), 2))));
				S0123 = 0.5 * fabs(sqrt((pow((D03[1] * D12[2] - D03[2] * D12[1]), 2)) + (pow((D03[2] * D12[0] - D03[0] * D12[2]), 2)) + (pow((D03[0] * D12[1] - D03[1] * D12[0]), 2))));
				S4567 = 0.5 * fabs(sqrt((pow((D47[1] * D56[2] - D47[2] * D56[1]), 2)) + (pow((D47[2] * D56[0] - D47[0] * D56[2]), 2)) + (pow((D47[0] * D56[1] - D47[1] * D56[0]), 2))));
				fprintf(stdout, "\n%e %e %e %e %e %e \n", S1357, S0145, S0246, S2367, S0123, S4567);

				/*
				 * NORMAL VECTOR definition Normal vector of
				 * a plane with known two vectors can be
				 * calculated by finding the cross product of
				 * these two vectors. Here, AC X BD will give
				 * the normal vector to the plane where AC
				 * and BD lie. Also, X[0] = i component X[1]
				 * = j component X[2] = k component
				 * 
				 * These vector are constructed in a way to have
				 * their norms equal to the surface S to
				 * which they are perpendicular e.g.  M = 0.5
				 * ( AC X BD )
				 */
				int             x = 0, y = 0;
				for (n = 0; n < 3; n++) {
					if (n == 0) {
						x = 1;
						y = 2;
					}
					if (n == 1) {
						x = 2;
						y = 0;
					}
					if (n == 2) {
						x = 0;
						y = 1;
					}
					N1357[n] = 0.5 * (D17[x] * D35[y] - D17[y] * D35[x]);
					N0145[n] = -0.5 * (D05[x] * D14[y] - D05[y] * D14[x]);
					N0246[n] = -0.5 * (D06[x] * D24[y] - D06[y] * D24[x]);
					N2367[n] = 0.5 * (D27[x] * D36[y] - D27[y] * D36[x]);
					N0123[n] = -0.5 * (D03[x] * D12[y] - D03[y] * D12[x]);
					N4567[n] = 0.5 * (D47[x] * D56[y] - D47[y] * D56[x]);
				}

				/*
				 * volume definition Ref:  Jeffrey Grandy,
				 * Efficient Computation of Volume of
				 * Hexahedral Cells, Lawrence Livermore
				 * National Laboratory
				 */

				for (n = 0; n < 3; n++) {
					M0167[n] = D17[n] + D06[n];
					M0257[n] = D05[n] + D27[n];
					M0347[n] = D03[n] + D47[n];
				}
				V1 = fabs(M0167[0] * (D27[1] * D03[2] - D27[2] * D03[1]) - D27[0] * (M0167[1] * D03[2] - M0167[2] * D03[1]) + D03[0] * (M0167[1] * D27[2] - M0167[2] * D27[1]));
				V2 = fabs(D06[0] * (M0257[1] * D47[2] - M0257[2] * D47[1]) - M0257[0] * (D06[1] * D47[2] - D06[2] * D47[1]) + D47[0] * (D06[1] * M0257[2] - D06[2] * M0257[1]));
				V3 = fabs(D17[0] * (D05[1] * M0347[2] - D05[2] * M0347[1]) - D05[0] * (D17[1] * M0347[2] - D17[2] * M0347[1]) + M0347[0] * (D17[1] * D05[2] - D17[2] * D05[1]));
				V = (V1 + V2 + V3) / 12;

				/* testing volume */
				fprintf(stdout, "%e \n", V);

				/*
			        nS_over_V[IDX(i,j,k)].xp.x = N1357[0]/V;
			        nS_over_V[IDX(i,j,k)].xp.y = N1357[1]/V;
			        nS_over_V[IDX(i,j,k)].xp.z = N1357[2]/V;
			       
			        nS_over_V[IDX(i,j,k)].xm.x = N0246[0]/V;
			        nS_over_V[IDX(i,j,k)].xm.y = N0246[1]/V;
			        nS_over_V[IDX(i,j,k)].xm.z = N0246[2]/V;
			       
			        nS_over_V[IDX(i,j,k)].yp.x = N2367[0]/V;
			        nS_over_V[IDX(i,j,k)].yp.y = N2367[1]/V;
			        nS_over_V[IDX(i,j,k)].yp.z = N2367[2]/V;
			       
			        nS_over_V[IDX(i,j,k)].ym.x = N0145[0]/V;
			        nS_over_V[IDX(i,j,k)].ym.y = N0145[1]/V;
			        nS_over_V[IDX(i,j,k)].ym.z = N0145[2]/V;
			       
			        nS_over_V[IDX(i,j,k)].zp.x = N4567[0]/V;
			        nS_over_V[IDX(i,j,k)].zp.y = N4567[1]/V;
			        nS_over_V[IDX(i,j,k)].zp.z = N4567[2]/V;
			       
			        nS_over_V[IDX(i,j,k)].zm.x = N0123[0]/V;
			        nS_over_V[IDX(i,j,k)].zm.y = N0123[1]/V;
			        nS_over_V[IDX(i,j,k)].zm.z = N0123[2]/V;
			        */

				center_V[IDX(i, j, k)].x = P8[0];
				center_V[IDX(i, j, k)].y = P8[1];
				center_V[IDX(i, j, k)].z = P8[2];


				for (pp = 0; pp < NPOP; pp++) {
					coeff_xp[IDX(i, j, k)].p[pp] = (N1357[0] * c[pp].x + N1357[1] * c[pp].y + N1357[2] * c[pp].z) / V;
					coeff_xm[IDX(i, j, k)].p[pp] = (N0246[0] * c[pp].x + N0246[1] * c[pp].y + N0246[2] * c[pp].z) / V;
					coeff_yp[IDX(i, j, k)].p[pp] = (N2367[0] * c[pp].x + N2367[1] * c[pp].y + N2367[2] * c[pp].z) / V;
					coeff_ym[IDX(i, j, k)].p[pp] = (N0145[0] * c[pp].x + N0145[1] * c[pp].y + N0145[2] * c[pp].z) / V;
					coeff_zp[IDX(i, j, k)].p[pp] = (N4567[0] * c[pp].x + N4567[1] * c[pp].y + N4567[2] * c[pp].z) / V;
					coeff_zm[IDX(i, j, k)].p[pp] = (N0123[0] * c[pp].x + N0123[1] * c[pp].y + N0123[2] * c[pp].z) / V;

				}

			}
		}
	}			/* for i , j , k */




}
