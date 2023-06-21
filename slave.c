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
__thread_local_share int shared_indx;
//csr_spmv + csr_precondition_spmv
__thread_local double matrix_1[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local double matrix_2[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local double vector[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local double result[dataBufferSize] __attribute__ ((aligned(64))); // gAPtr

//pcg_precondition_csr
__thread_local_share double wAPtr[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local_share double result_air[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local double preD[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local double rAPtr[dataBufferSize] __attribute__ ((aligned(64)));
__thread_local double wAPtr_local[dataBufferSize] __attribute__ ((aligned(64)));

void slave_csr_spmv(task_csr_spmv * tp)
{
	task_csr_spmv slave_task;
	double * matrix = nullptr;
	CRTS_dma_iget(&slave_task, tp, sizeof(task_csr_spmv), &DMARply);
	DMARplyCount++;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);
	// CRTS_tid
	
	//csr 信息
	int rows = slave_task.rows;
	int *row_off = slave_task.row_off
	int *cols = slave_task.cols;
	double *data = slave_task.data;
	int data_size = slave_task.data_size;

	//vector 信息
	int length = slave_task.length;
	double * vec = slave_task.vec;

	//全局信息
	int srow = slave_task.srow;
	int max_entry = slave_task.max_entry;

	//持久的信息
	athread_memcpy_sldm(&shared_indx,&tp->cur_indx,sizeof(int),MEM_TO_LDM);
	CRTS_dma_iget(&vector,slave_task.vec,length * sizeof(double), &DMARply);
	//双缓存矩阵

	double * start_addr = data + (row_off[CRTS_tid])*sizeof(double);
	int start_indx = CRTS_tid;
	int end_indx = min(start_indx+srow,rows);
	int gap = end_indx - start_indx;
	int length_iget = (row_off[end_indx] - row_off[start_indx])*sizeof(double);
	CRTS_dma_iget(&matrix_2,slave_task.data,length_iget, &DMARply);
	DMARplyCount += 2;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);
	int flag = 2;

	while (shared_indx <rows)
	{
		

		//计算稀疏矩阵
		if (flag ==2)
		{
			matrix = matrix_2;
			//确定读取空间
			CRTS_smutex_lock_array();
			start_addr = data + (row_off[shared_indx])*sizeof(double)；
			start_indx = shared_indx;
			end_indx = min(start_indx+srow,rows);
			shared_indx = end_indx;
			CRTS_smutex_unlock_array();

			// 读取部分稀疏矩阵
			gap = end_indx - start_indx;
			length_iget = (row_off[end_indx] - row_off[start_indx])*sizeof(double);
			CRTS_dma_iget(&matrix_1,start_addr,length_iget,&DMARply)
			DMARplyCount += 1;
			CRTS_dma_wait_value(&DMARply, DMARplyCount);
			//更改标志
			flag = 1;
		}
		else
		{
			matrix = matrix_1;
			//确定读取空间
			CRTS_smutex_lock_array();
			start_addr = data + (row_off[shared_indx])*sizeof(double);
			start_indx = shared_indx;
			end_indx = min(start_indx+srow,rows);
			shared_indx = end_indx;
			CRTS_smutex_unlock_array();

			// 读取部分稀疏矩阵
			gap = end_indx - start_indx;
			length_iget = (row_off[end_indx] - row_off[start_indx])*sizeof(double);
			CRTS_dma_iget(&matrix_2,start_addr,length_iget,&DMARply)
			DMARplyCount += 1;
			CRTS_dma_wait_value(&DMARply, DMARplyCount);
			//更改标志
			flag = 2;
		}
		for (int i=0;i<gap;i++)
		{
			int start = row_off[start_indx+i];
			int num = row_off[start_indx+i+1] - row_off[start_indx+i];
			double tmp = 0.0
			for (int j=0;j<num;j++)
			{
				tmp+= vector[cols[start+j]]*matrix[cols[start+j]];
			}
			result[start_indx+i] = tmp;
		}
	}
}


void slave_csr_precondition_spmv(task_pcg_precondition_csr * tp)
{
	task_pcg_precondition_csr slave_task;
	CRTS_dma_iget(&slave_task, tp, sizeof(task_pcg_precondition_csr), &DMARply);
	DMARplyCount++;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);

	int length_preD = (slave_task.length_preD +63)/64;
	int length_rAPtr = length_preD;

	preD_start_addr = CRTS_tid*length_preD + slave_task.preD;
	rAPtr_start_addr = CRTS_tid*length_rAPtr + slave_task.rAPtr;

	length = min(length_preD,slave_task.length_preD-CRTS_tid*length_preD);

	CRTS_dma_iget(&preD,slave_task.preD,length * sizeof(double), &DMARply);
	CRTS_dma_iget(&rAPtr,slave_task.rAPtr,length * sizeof(double), &DMARply);
	
	//simd 
	doublev8 va,vb,vc;
	int i =0;
	for (;i<(length/8)*8;i+=8)
	{
		simd_load(va,preD+i);
		simd_load(vb,rAPtr+i);
		vc = va * vb;
		simd_store(vc,wAPtr_local+i);
	}
	for (;i<length;i++)
	{
		wAPtr_local[i] = preD[i]*rAPtr[i];
	}
	CRTS_sldm_get(wAPtr+length_preD*CRTS_tid,wAPtr_local+length_preD*CRTS_tid,length);



	double * val = nullptr;
	
	// CRTS_tid
	
	//csr 信息
	int rows = slave_task.rows;
	int *row_off = slave_task.row_off
	int *cols = slave_task.cols;

	//vector1 信息
	int length = slave_task.length_preD;
	double * vec = slave_task.vec;

	int size = slave_task.size;
	double *val = slave_task.val

	//全局信息
	int srow = slave_task.srow;
	int max_entry = slave_task.max_entry;

	//持久的信息
	athread_memcpy_sldm(&shared_indx,&tp->cur_indx,sizeof(int),MEM_TO_LDM);
	CRTS_sldm_get(&vector,wAPtr,length * sizeof(double));
	//双缓存矩阵

	double * start_addr = val + (row_off[CRTS_tid])*sizeof(double);
	int start_indx = CRTS_tid;
	int end_indx = min(start_indx+srow,rows);
	int gap = end_indx - start_indx;
	int length_iget = (row_off[end_indx] - row_off[start_indx])*sizeof(double);
	CRTS_dma_iget(&matrix_2,slave_task.val,length_iget, &DMARply);
	DMARplyCount += 2;
	CRTS_dma_wait_value(&DMARply, DMARplyCount);
	int flag = 2;

	while (shared_indx <rows)
	{
		//计算稀疏矩阵
		if (flag ==2)
		{
			matrix = matrix_2;
			//确定读取空间
			CRTS_smutex_lock_array();
			start_addr = val + (row_off[shared_indx])*sizeof(double)；
			start_indx = shared_indx;
			end_indx = min(start_indx+srow,rows);
			shared_indx = end_indx;
			CRTS_smutex_unlock_array();

			// 读取部分稀疏矩阵
			gap = end_indx - start_indx;
			length_iget = (row_off[end_indx] - row_off[start_indx])*sizeof(double);
			CRTS_dma_iget(&matrix_1,start_addr,length_iget,&DMARply)
			DMARplyCount += 1;
			CRTS_dma_wait_value(&DMARply, DMARplyCount);
			//更改标志
			flag = 1;
		}
		else
		{
			matrix = matrix_1;
			//确定读取空间
			CRTS_smutex_lock_array();
			start_addr = val + (row_off[shared_indx])*sizeof(double);
			start_indx = shared_indx;
			end_indx = min(start_indx+srow,rows);
			shared_indx = end_indx;
			CRTS_smutex_unlock_array();

			// 读取部分稀疏矩阵
			gap = end_indx - start_indx;
			length_iget = (row_off[end_indx] - row_off[start_indx])*sizeof(double);
			CRTS_dma_iget(&matrix_2,start_addr,length_iget,&DMARply)
			DMARplyCount += 1;
			CRTS_dma_wait_value(&DMARply, DMARplyCount);
			//更改标志
			flag = 2;
		}
		for (int i=0;i<gap;i++)
		{
			int start = row_off[start_indx+i];
			int num = row_off[start_indx+i+1] - row_off[start_indx+i];
			double tmp = 0.0
			for (int j=0;j<num;j++)
			{
				tmp+= vector[cols[start+j]]*matrix[cols[start+j]];
			}
			result[start_indx+i] = tmp;
			result_air[start_indx+i] = tmp;
			CRTS_sldm_get(result_air+start_indx+i,result+start_indx+i,sizeof(double));
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
