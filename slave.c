#include <slave.h>
#include "pcg_def.h"

#include <crts.h>

typedef struct{
	double *p;
	double *z;
	double beta;
	int cells;
} Para;

#define dataBufferSize 2000
__thread_local crts_rply_t DMARply = 0;
__thread_local unsigned int DMARplyCount = 0;
// __thread_local double p[dataBufferSize] __attribute__ ((aligned(64)));
// __thread_local double z[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local double sw_mat[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local double x[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local double y[dataBufferSize] __attribute__ ((aligned(64)));
void slave_spmv(task_pool * tp)
{
	task_pool slave_tp;
	CRTS_dma_iget(&slave_tp, tp, sizeof(task_pool), &DMARply);
	DMARplyCount++;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);
	// CRTS_tid

	int srow = slave_tp._srow;
	int indx = slave_tp._indx;
	int rows = slave_tp.rows;
	int *cols = slave_tp.cols;

	int * start = slave_tp.row_off[indx];
	int * end = slave_tp.row_off[min(rows,indx+srow)]

	int length = slave_tp.length
	double * vec = slave_tp.vec_data

	CRTS_dma_iget(&sw_mat, start, (end-start) * sizeof(double), &DMARply);
	CRTS_dma_iget(&x, vec, length * sizeof(double), &DMARply);
	DMARplyCount += 2;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);

	for (int i=0;i<min(rows,indx+srow) - indx;i++)
	{
		int start_ = start[i] 
		num = slave_tp.row_off[indx+i+1] - slave_tp.row_off[indx+i]
		double temp = 0;
		for (int j=0;j<num;j++)
		{
			temp += x[cols[j+start_]]*sw_mat[j+start_]
		}
		y[i] = temp
	}



}

	






}
void slave_example(Para* para){
	Para slavePara;
	//接收结构体数据
	CRTS_dma_iget(&slavePara, para, sizeof(Para), &DMARply);
	DMARplyCount++;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);
	double beta = slavePara.beta;
	int cells = slavePara.cells;
	
	//计算从核接收数组数据长度和接收位置
	int len = cells / 64;
	int rest = cells % 64;
	int addr;
	if(CRTS_tid < rest){
		len++;
		addr = CRTS_tid * len;
	}else{
		addr = CRTS_tid * len + rest;
	}
	//接收数组数据
	CRTS_dma_iget(&p, slavePara.p + addr, len * sizeof(double), &DMARply);
	CRTS_dma_iget(&z, slavePara.z + addr, len * sizeof(double), &DMARply);
	DMARplyCount += 2;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);
			
	//计算
	int i = 0;
	for(; i < len; i++){
		p[i] = z[i] + beta * p[i];
	}
	//传回计算结果
	CRTS_dma_iput(slavePara.p+addr, &p, len * sizeof(double), &DMARply);
	DMARplyCount++;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);
}
