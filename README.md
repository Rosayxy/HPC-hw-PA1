## 重要
1. 性能分
   1. 前两组单线程调用std::sort 稳过
   2. 第三组用 **std::sort** `srun -n 16 odd_even_sort 10000 data/10000.dat` 算平均值应该不会超
   3. 第四组换 radix sort `srun -n 20 odd_even_sort 100000 data/100000.dat` 大部分在1.0几毫秒，平均不会超
   4. 第五组 radix sort `srun -n 56 --cpu_bind=sockets odd_even_sort 1000000 data/1000000.dat` radix sort 在5.5ms 以内
   5. 第六组`srun -N 2 -n 56 odd_even_sort 10000000 data/10000000.dat` 50 毫秒以下 -n 50 62 毫秒
   6. 第七组 `srun -N 2 -n 56 odd_even_sort 100000000 data/100000000.dat` 650毫秒之内 -n 56 700-740 