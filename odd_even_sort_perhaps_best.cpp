#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mpi.h>
#include <vector>
#include "worker.h"
//# define debug 1
// todo 判断什么时候连续两次都没有交换 不用cnt 了 改成 prev_flag 吧
// todo 不行的话试一下 BSend 吧
#define array(a, b, c, block_len) (((size_t)(c) >= (block_len)) ? ((b)[c - block_len]) : ((a)[c]))
void Worker::sort()
{
  /** Your code ... */
  // you can use variables in class Worker: n, nprocs, rank, block_len, data
  // 拿取自己的块并且排序
//std::sort(data, data + block_len);
  // bucket sort
    # ifdef debug
  char c=rank+'0';
  char filename[]={c,'.','t','x','t','\0'};
  FILE* f=fopen(filename,"w");
  fprintf(f,"at sort start: %d ",rank);
  fflush(f);
  for(int i=0;(size_t)i<block_len;i++){
    fprintf(f,"%f ",data[i]);
  }
  fprintf(f,"\n");
  fflush(f);
  # endif
  // todo 测试负数上的兼容性
    std::vector<float> bucket[256];
  for(size_t i=0;i<block_len;i++){
    int *p=(int*)&data[i];
    *p=((*p)>>31 & 0x1)? ~(*p) : (*p) | 0x80000000;
  }
    for (int i=0; i<4; i++) {
        for (int j=0; (size_t)j<block_len; j++)
            bucket[(*(int*)(&data[j]))>>(i*8) & 0xff].push_back(data[j]);
        int count = 0;
        for (int j=0; j<256; j++) {
            for (int k=0; k<(int)bucket[j].size(); k++)
                data[count++] = bucket[j][k];
            bucket[j].clear();
        }
    }
    for(size_t i=0;i<block_len;i++){
        int *p=(int*)&data[i];
        *p=((*p)>>31 & 0x1)? (*p) & 0x7fffffff : ~(*p);
    }
  if (nprocs == 1)
  {
    return;
  }
  // int cnt = 0;  // previously is sorted
  // int flag = 0;
  // 复用worker.h 中的 last_rank, out_of_range，分别表示本组是否需要交换和整体是否需要交换
  // 互相交换元素，如果是偶数块，发现需要排序就像前一个奇数块发送flag,像后一个奇数块发自己所有块，
  //
  float buf = 0;
  float *tmp = new float[block_len+1]; // todo 为了兼容性改一下
  float* buf1=new float[block_len+1];
  //MPI_Request request;
  //MPI_Request request1;
  for (int ii = 0; ii < (nprocs ) / 2; ii++)
  {
    // 偶数排序
   // MPI_Barrier(MPI_COMM_WORLD);
    if (rank % 2 == 0)
    {
      // int MPI_Sendrecv(MPICH2_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype,int dest, int sendtag,void *recvbuf, int recvcount, MPI_Datatype recvtype,int source, int recvtag,MPI_Comm comm, MPI_Status *status)
      MPI_Sendrecv(data + block_len - 1, 1, MPI_FLOAT, rank + 1, 0, &buf, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (data[block_len - 1] > buf)
      {
        // 接收奇数块的数据
        MPI_Sendrecv(data,block_len,MPI_FLOAT,rank+1,0,buf1, block_len, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // merge
        //std::merge(data,data+block_len,buf1,buf1+block_len,tmp);
        float *start_ptr=data,*start1_ptr=buf1,*res=tmp,*last=data+block_len,*last1=buf1+block_len;
        # ifdef debug
        fprintf(f,"at even merge start: %d ",rank);
        fprintf(f,"%p %p %p %p %p %p %f %f\n",(void*)start_ptr,(void*)start1_ptr,(void*)res,(void*)last,(void*)last1,(void*)tmp,*start_ptr,*start1_ptr);
        fflush(f);
        # endif
        size_t cnt=0;
        // 只用排序前一半
        while((cnt<block_len)&&start_ptr<last&&start1_ptr<last1){
            if(*start_ptr<*start1_ptr){
              # ifdef debug
              fprintf(f,"inside while...\n");
              fprintf(f,"%p %p %f %f \n",(void*)start_ptr,(void*)res,*start_ptr,*res);
              fflush(f);
              # endif
                *res=*start_ptr;
                start_ptr++;
            }else{
              # ifdef debug
              fprintf(f,"inside while...\n");
              fprintf(f,"%p %p %f %f \n",(void*)start1_ptr,(void*)res,*start1_ptr,*res);
              fflush(f);
              # endif
                *res=*start1_ptr;
                start1_ptr++;
            }
            res++;
            cnt++;
        }
        // 拷贝回原先data
       // memcpy(data,tmp,block_len*sizeof(float));
        # ifdef debug
        fprintf(f,"at even merge end: %d ",rank);
        for(int i=0;(size_t)i<block_len;i++){
            fprintf(f,"%f ",tmp[i]);
        }
        fprintf(f,"\n");
        fflush(f);
        # endif
        float* tmmp=data;
        data=tmp;
        tmp=tmmp;
        # ifdef debug
        fprintf(f,"at even swap end: %d ",rank);
        fprintf(f,"%p %p\n",(void*)data,(void*)tmp);
        # endif
      }
    }
    else
    {
      MPI_Sendrecv(data, 1, MPI_FLOAT, rank - 1, 0, &buf, 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (data[0] < buf)
      {
        MPI_Sendrecv(data,block_len,MPI_FLOAT,rank-1,0,buf1, block_len, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // 反向 merge
        float *start=data+block_len-1,*start1=buf1+block_len-1,*res=tmp+block_len-1,*last=data-1,*last1=buf1-1;
        size_t cnt=0;
        # ifdef debug
        fprintf(f,"at merge start: %d ",rank);
        fprintf(f,"%p %p %p %p %p %p\n",(void*)start,(void*)start1,(void*)res,(void*)last,(void*)last1,(void*)tmp);
        fflush(f);
        # endif
        while((cnt<block_len)&&(start>last)&&start1>last1){
            if(*start>*start1){
                *res=*start;
                start--;
            }else{
                *res=*start1;
                start1--;
            }
            res--;
            cnt++;
        }
        # ifdef debug
        fprintf(f,"at even merge end: %d ",rank);
        for(int i=0;(size_t)i<block_len;i++){
            fprintf(f,"%f ",tmp[i]);
        }
        fprintf(f,"\n");
        fflush(f);
        # endif
        float* tmmp=data;
        data=tmp;
        tmp=tmmp;
        # ifdef debug
        fprintf(f,"at even swap end: %d ",rank);
        fprintf(f,"%p %p\n",(void*)data,(void*)tmp);
        # endif
      }
    }
   // MPI_Barrier(MPI_COMM_WORLD);
    // 奇数排序
    if (rank % 2 == 1)
    {
      // 和自己后一个偶数块绑定
      if (rank != nprocs - 1)
      {
        MPI_Sendrecv(data + block_len - 1, 1, MPI_FLOAT, rank + 1, 0, &buf, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (data[block_len - 1] > buf)
        {
          // merge
          // 接收奇数块的数据
        MPI_Sendrecv(data,block_len,MPI_FLOAT,rank+1,0,buf1, block_len, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // merge
        //std::merge(data,data+block_len,buf1,buf1+block_len,tmp);
        float *start_ptr=data,*start1_ptr=buf1,*res=tmp,*last=data+block_len,*last1=buf1+block_len;
        size_t cnt=0;
        # ifdef debug
        fprintf(f,"at odd merge start: %d ",rank);
        fprintf(f,"%p %p %p %p %p %p\n",(void*)start_ptr,(void*)start1_ptr,(void*)res,(void*)last,(void*)last1,(void*)tmp);
        fflush(f);
        # endif
        // 只用排序前一半
        while((cnt<block_len)&&start_ptr<last&&start1_ptr<last1){
            if(*start_ptr<*start1_ptr){
                *res=*start_ptr;
                start_ptr++;
            }else{
                *res=*start1_ptr;
                start1_ptr++;
            }
            res++;
            cnt++;
        }
         # ifdef debug
        fprintf(f,"at odd merge end: %d ",rank);
        for(int i=0;(size_t)i<block_len;i++){
            fprintf(f,"%f ",tmp[i]);
        }
        fprintf(f,"\n");
        fflush(f);
        # endif
        // 拷贝回原先data
       // memcpy(data,tmp,block_len*sizeof(float));
        float* tmmp=data;
        data=tmp;
        tmp=tmmp;
        # ifdef debug
        fprintf(f,"at odd swap end: %d ",rank);
        fprintf(f,"%p %p\n",(void*)data,(void*)tmp);
        # endif
        }
      }
    }
    else if (rank % 2 == 0)
    {
      if (rank != 0)
      {
        // 如果需要交换 则像自己前一个奇数块发送自己所有数据
        MPI_Sendrecv(data, 1, MPI_FLOAT, rank - 1, 0, &buf, 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (data[0] < buf)
        {
        MPI_Sendrecv(data,block_len,MPI_FLOAT,rank-1,0,buf1, block_len, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // 反向 merge
        float *start=data+block_len-1,*start1=buf1+block_len-1,*res=tmp+block_len-1,*last=data-1,*last1=buf1-1;
        size_t cnt=0;
        # ifdef debug
        fprintf(f,"at odd merge start: %d ",rank);
        fprintf(f,"%p %p %p %p %p %p\n",(void*)start,(void*)start1,(void*)res,(void*)last,(void*)last1,(void*)tmp);
        fflush(f);
        # endif
        while((cnt<block_len)&&(start>last)&&start1>last1){
            if(*start>*start1){
                *res=*start;
                start--;
            }else{
                *res=*start1;
                start1--;
            }
            res--;
            cnt++;
        }
        # ifdef debug
        fprintf(f,"at odd merge end: %d ",rank);
        for(int i=0;(size_t)i<block_len;i++){
            fprintf(f,"%f ",tmp[i]);
        }
        fprintf(f,"\n");
        fflush(f);
        # endif
        float* tmmp=data;
        data=tmp;
        tmp=tmmp;
        # ifdef debug
        fprintf(f,"at odd swap end: %d ",rank);
        fprintf(f,"%p %p\n",(void*)data,(void*)tmp);
        # endif
        }
      }
    }
  }

//  MPI_Barrier(MPI_COMM_WORLD);
    if (rank % 2 == 0)
    {
      // int MPI_Sendrecv(MPICH2_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype,int dest, int sendtag,void *recvbuf, int recvcount, MPI_Datatype recvtype,int source, int recvtag,MPI_Comm comm, MPI_Status *status)
      MPI_Sendrecv(data + block_len - 1, 1, MPI_FLOAT, rank + 1, 0, &buf, 1, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (data[block_len - 1] > buf)
      {
        // 接收奇数块的数据
        MPI_Sendrecv(data,block_len,MPI_FLOAT,rank+1,0,buf1, block_len, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // merge
        //std::merge(data,data+block_len,buf1,buf1+block_len,tmp);
        float *start_ptr=data,*start1_ptr=buf1,*res=tmp,*last=data+block_len,*last1=buf1+block_len;
        size_t cnt=0;
        // 只用排序前一半
        while((cnt<block_len)&&start_ptr<last&&start1_ptr<last1){
            if(*start_ptr<*start1_ptr){
                *res=*start_ptr;
                start_ptr++;
            }else{
                *res=*start1_ptr;
                start1_ptr++;
            }
            res++;
            cnt++;
        }
        // 拷贝回原先data
       // memcpy(data,tmp,block_len*sizeof(float));
        float* tmmp=data;
        data=tmp;
        tmp=tmmp;
      }
    }
    else
    {
      MPI_Sendrecv(data, 1, MPI_FLOAT, rank - 1, 0, &buf, 1, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if (data[0] < buf)
      {
        MPI_Sendrecv(data,block_len,MPI_FLOAT,rank-1,0,buf1, block_len, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // 反向 merge
        float *start=data+block_len-1,*start1=buf1+block_len-1,*res=tmp+block_len-1,*last=data-1,*last1=buf1-1;
        size_t cnt=0;
        while((cnt<block_len)&&(start>last)&&start1>last1){
            if(*start>*start1){
                *res=*start;
                start--;
            }else{
                *res=*start1;
                start1--;
            }
            res--;
            cnt++;
        }float* tmmp=data;
        data=tmp;
        tmp=tmmp;
      }
    }
  //  MPI_Barrier(MPI_COMM_WORLD);
}
