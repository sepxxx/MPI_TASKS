
#include <mpi.h>
//#include <iostream>
#include <random>
//#include<cstdlib> 
using namespace std;
void t1() {
	//MPI_COMM_WORLD 
	//int rank, size;
	//MPI_Comm_rank();
	MPI_Init(NULL, NULL);
	printf("Hello world!\n");
	MPI_Finalize();
}

void t2() {
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const int n = 10;
	int* data = new int[n];
	MPI_Status st;
	if (rank == 0) {
		for (int i = 0; i < n; i++)
			data[i] = i*i;
	}
	MPI_Bcast(data, n, MPI_INT, 0, MPI_COMM_WORLD);

	printf("\nRANK: %d\n", rank);
	for (int i = 0; i < n; i++)
		printf(" %d ", data[i]);

	int col_each_proc, ost = 0;
	col_each_proc = n / size;
	if (n % size)
		ost = n - col_each_proc * size;
	//printf("%d %d\n", col_each_proc, ost);
	
	int istart = rank * col_each_proc;
	int procmax = data[istart];
	int globmax;
	//j сколько раз крутимся
	//i для обращений по индексам
	for (int i = istart, j = 0; j < col_each_proc; j++, i++) {
		//printf(" %d ", i);
		if (data[i] > procmax)
			procmax = data[i];
	}
	//процессу с таким индексом отдаем остаток
	if ((rank == size - 1) && ost) {
		//printf("OST: %d\n", ost);
		for (int i = (rank + 1) * col_each_proc, j = 0; j < ost; j++, i++) {
			//printf(" %d ", i);
			if (data[i] > procmax)
				procmax = data[i];
		}
	}

	MPI_Reduce(&procmax, &globmax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	if(rank==0)
		printf("\nMAXIMUM: %d\n", globmax);
	MPI_Finalize();
		
}

void t3() {
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int Npoints = 1e6;
	int col_each_proc, ost = 0;
	col_each_proc = Npoints / size;
	if (Npoints % size)
		ost = Npoints  - col_each_proc * size;

	int Nincircle = 0, Ncom = 0;
	double lbound = -1;
	double ubound = 1;
	uniform_real_distribution<double> unif(lbound, ubound);
	default_random_engine re;
	int end = col_each_proc;
	if (ost && rank == size - 1)
		end += ost;
		
	for (int j = 0; j < end; j++) {
		double x = unif(re);
		double y = unif(re);
		if((x * x + y * y) < 1)
			Nincircle++;
	}


	MPI_Reduce(&Nincircle, &Ncom, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		double pi = Ncom * 4.0 / Npoints;
		printf("RES: %lf\n", pi);
	}
	MPI_Finalize();

}

void t4() {
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double lbound = -1;
	double ubound = 1;
	uniform_real_distribution<double> unif(lbound, ubound);
	default_random_engine re;
	const int n = 1200;
	int locn = n / size;
	double locdata[n];
	double data[n];

	if (rank == 0) {
		for (int j = 0; j < n; j++) {
			data[j] = unif(re);
		}
	}
	const int sendcounts[4] = { locn, locn, locn, locn };
	const int displs[4] = { 0, locn, 2*locn, 3*locn };
	//MPI_Scatter(data, locn, MPI_DOUBLE, &locdata[0], locn, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(data, sendcounts, displs, MPI_DOUBLE, locdata, locn, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		double locsum = 0;
		int loc_plus = 0;
		for (int i = 0; i < locn; i++) {
			if (locdata[i] > 0) {
				locsum += locdata[i];
				loc_plus++;
			}

		}
		//printf("LOCPLUS: %d\n", loc_plus);
		int rcol;
		double rsum;
		MPI_Reduce(&loc_plus, &rcol, 1,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&locsum, &rsum, 1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (rank == 0) {
			printf("rcol: %d\n", rcol);
			printf("rsum: %lf\n", rsum);
			printf("mean: %lf\n", rsum/rcol);
		}

	MPI_Finalize();

}

void t5() {
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double lbound = 0;
	double ubound = 10;
	uniform_real_distribution<double> unif(lbound, ubound);
	default_random_engine re;
	const int n = 1000;
	int locn = n / size;
	double ld1[n];
	double ld2[n];

	double d1[n];
	double d2[n];

	if (rank == 0) {
		for (int j = 0; j < n; j++) {
			d1[j] = unif(re);
			d2[j] = unif(re);
		}
	}
	const int sendcounts[4] = { locn, locn, locn, locn };
	const int displs[4] = { 0, locn, 2 * locn, 3 * locn };
	MPI_Scatterv(d1, sendcounts, displs, MPI_DOUBLE, ld1, locn, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(d2, sendcounts, displs, MPI_DOUBLE, ld2, locn, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//MPI_Scatter(&d1[0], locn, MPI_DOUBLE, &ld1[0], locn, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Scatter(&d2[0], locn, MPI_DOUBLE, &ld2[0], locn, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double locsum = 0;
	for (int i = 0; i < locn; i++) {
		locsum += ld1[i] * ld2[i];
	}
	double rsum;
	MPI_Reduce(&locsum, &rsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		printf("rsum: %lf\n", rsum);
	}

	MPI_Finalize();

}

void t6() {
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double lbound = 0;
	double ubound = 10;
	
	uniform_real_distribution<double> unif(lbound, ubound);
	default_random_engine re;
	const int n = 12;
	//glob
	double** pa = (double**)malloc(n * sizeof(double*));
	double* va = (double*)malloc(n * n * sizeof(double));
	for (int i = 0; i < n; i++)
		pa[i] = va + n * i;
	//local
	double** pa_l = (double**)malloc(n * sizeof(double*));
	double* va_l = (double*)malloc(n * n * sizeof(double));
	for (int i = 0; i < n; i++)
		pa_l[i] = va_l + n * i;

	//printf("RANK %d\n", rank);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			pa[i][j] = unif(re);
			//printf(" %lf ", pa[i][j]);
		}
		//printf("\n");
	}
	
	int col_rows = n / size;
	//printf("NROWS %d \n", n * col_rows);

	MPI_Scatter(&pa[0][0], n * col_rows, MPI_DOUBLE, &pa_l[0][0], n * col_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	

	//printf("RANK %d\n", rank);
	double maxmin=-9999, min= pa_l[0][0];
	for (int i = 0; i < col_rows; i++) {
		for (int j = 0; j < n; j++) {
			if (pa_l[i][j] < min)
				min = pa_l[i][j];
		}
		if (maxmin < min)
			maxmin = min;
		//printf("\n");
	}
	double globminmax;
	MPI_Reduce(&maxmin, &globminmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		printf("SEDLOVAYA: %lf\n", globminmax);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				printf(" %.1lf ", pa[i][j]);
			}
			printf("\n");
		}
	}
	MPI_Finalize();

}

void t8() {
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const int n = 1e2;
	int ln = n / size;
	int* ldata = new int[ln];
	MPI_Status st;
	if (rank == 0) {
		int* data = new int[n];
		for (int i = 0; i < n;i++)
			data[i] = i;

		printf("\n-----------------------------ARRAY FOR SCATTER---------------------------------\n");
		for (int i = 0; i < n; i++)
			printf(" %d ", data[i]);
		
		for (int i = 1; i < size; i++) {
			MPI_Send(data + i * ln, ln, MPI_INT, i, 666, MPI_COMM_WORLD);
			//printf(" %d ", *(data + i *ln);
		}
		for (int i = 0; i < ln;i++)
			ldata[i] = data[i];
	}
	if(rank)
		MPI_Recv(ldata, ln, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
	
	printf("\nRANK: %d\n", rank);
	for (int i = 0; i < ln; i++)
		printf(" %d ", ldata[i]);

	if (rank!=3) {
		for (int i = 0; i < size - 1; i++) {
			MPI_Send(ldata + i * ln, ln, MPI_INT, 3, 666, MPI_COMM_WORLD);
		}
	}
	if (rank == 3) {
		int* rdata = new int[n];
		for (int i = 0; i < size - 1; i++)
			MPI_Recv(rdata + i*ln, ln, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
		for (int i = 3 * ln,j=0; i < n, j<ln; i++, j++)
			rdata[i] = ldata[j];
		printf("\n-----------------------------RESULT OF GATHER---------------------------------\n");
		for (int i = 0; i < n; i++)
			printf(" %d ", rdata[i]);
	}
	MPI_Finalize();
}

void t9() {
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const int n = 36;
	srand(time(0));
	int data[n];
	int locn = n / size;
	int ldata[n];
	int resdata[n];

	for (int i = 0; i < n; i++) {
		//data[i] = rand() % 10;
		data[i] = i;

	}
	if (rank == 0) {
		for (int i = 0; i < n; i++)
		printf(" %d ", data[i]);
		printf("\n---------------------------------------------\n");
	}
	MPI_Scatter(data, locn, MPI_INT, ldata, locn, MPI_INT, 0, MPI_COMM_WORLD);
	
	//printf("RANK %d\n", rank);
	/*for (int i = 0; i < locn; i++) {
		printf(" %d ", ldata[i]);
	}*/
	for (int low = 0, high = locn - 1; low < high; low++, high--) {
		int tmp = ldata[low];
		ldata[low] = ldata[high];
		ldata[high] = tmp;
	}
	//printf("\n---------------------\n");
	//for (int i = 0; i < locn; i++) {
	//	printf(" %d ", ldata[i]);
	//}
	const int recvcounts[4] = { locn, locn, locn, locn };
	int ms = locn;
	const int displs[4] = { 3*ms, 2*ms, ms, 0 };
	MPI_Gatherv(ldata, locn, MPI_INT, resdata, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		for (int i = 0; i < n; i++) {
			printf(" %d ", resdata[i]);
		}
	}
	MPI_Finalize();
}

void t10() {
	MPI_Init(NULL, NULL);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	const int n = 120000;
	int* data = new int[n];
	int* ldata = new int[n];
	for (int i = 0; i < n; i++)
		data[i] = i;


	double starttime, endtime;
	MPI_Status st;
	//ДЛЯ BSEND
	int* buffer = new int[n+ MPI_BSEND_OVERHEAD];
	int bsize = sizeof(int) * (n+ MPI_BSEND_OVERHEAD);
	MPI_Buffer_attach(buffer, bsize);

	//
	if (rank == 0) {
		starttime = MPI_Wtime();
		//SEND 0.000941
		//MPI_Send(&data[0], n, MPI_INT, 1, 666, MPI_COMM_WORLD);
		//MPI_Recv(data, n, MPI_INT, 1, 777, MPI_COMM_WORLD, &st);
		//BSEND 0.001137
		//MPI_Bsend(&data[0], n, MPI_INT, 1, 666, MPI_COMM_WORLD);
		//MPI_Recv(data, n, MPI_INT, 1, 777, MPI_COMM_WORLD, &st);
		//RSEND 0.000995
		//начинается, если уже зарегистрирован соответствующий прием
		/*MPI_Rsend(&data[0], n, MPI_INT, 1, 666, MPI_COMM_WORLD);
		MPI_Recv(data, n, MPI_INT, 1, 777, MPI_COMM_WORLD, &st);*/
		//SSEND 0.001208
		//Адресат посылает источнику "квитанцию" - уведомление о завершении приема
		MPI_Ssend(&data[0], n, MPI_INT, 1, 666, MPI_COMM_WORLD);
		MPI_Recv(data, n, MPI_INT, 1, 777, MPI_COMM_WORLD, &st);

		endtime = MPI_Wtime();
	}
	if (rank == 1) {
		//MPI_Recv(&ldata[0], n, MPI_INT, 0, 666, MPI_COMM_WORLD, &st);
		//MPI_Send(ldata, n, MPI_INT, 0, 777, MPI_COMM_WORLD);

		//MPI_Recv(&ldata[0], n, MPI_INT, 0, 666, MPI_COMM_WORLD, &st);
		//MPI_Bsend(ldata, n, MPI_INT, 0, 777, MPI_COMM_WORLD);
		//MPI_Buffer_detach(&buffer, &bsize);

		//MPI_Recv(&ldata[0], n, MPI_INT, 0, 666, MPI_COMM_WORLD, &st);
		//MPI_Ssend(ldata, n, MPI_INT, 0, 777, MPI_COMM_WORLD);

		//MPI_Recv(&ldata[0], n, MPI_INT, 0, 666, MPI_COMM_WORLD, &st);
		//MPI_Rsend(ldata, n, MPI_INT, 0, 777, MPI_COMM_WORLD);

		MPI_Recv(&ldata[0], n, MPI_INT, 0, 666, MPI_COMM_WORLD, &st);
		MPI_Ssend(ldata, n, MPI_INT, 0, 777, MPI_COMM_WORLD);


	}
	if (rank == 0) {
		printf("That took %f seconds\n", endtime - starttime);
	}
	MPI_Finalize();

}


int main(int argc, char** argv) {
	//t1();
	//t2();
	//t3();
	//t4();
	//t5();
	//t6();
	//t8();
	//t9();
	//t10();


	///4 9 scatter v нормально распределить
	return 0;
}