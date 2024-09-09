// 在前一半在 data 后一半在 tmp 的情况下进行 merge
# include <cstdio>
# define block_len 5 

# define array(a,b,c) (((c)>=block_len)?((b)[c-block_len]):((a)[c]))
int main(){
float data[5]={2.3,3.6,7,9.1,10};
float tmp[5]={1.2,4.5,5.6,6.7,8.9};

if (data[block_len-1]>tmp[0]){
    //start merge
    int start=0,mid=block_len-1,end=2*block_len-1,start2=block_len;
    while(start<=mid&&start2<=end){
        if(array(data,tmp,start)<array(data,tmp,start2)){
            start++;
        }else{
            int val=array(data,tmp,start2);
            for(int i=start2;i>start;i--){
                array(data,tmp,i)=array(data,tmp,i-1);
            }
            array(data,tmp,start)=val;
            start++;
            start2++;
            mid++;
        }
    }
}
for(int i=0;i<block_len;i++){
    printf("%f ",data[i]);
}
for(int i=0;i<block_len;i++){
    printf("%f ",tmp[i]);
}
return 0;
}