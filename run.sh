#!/bin/bash
# 远程已经改成全 le 了
# run on 1 machine * 28 process, feel free to change it!
if [ "$2" -le 1000 ]; then
  srun -N 1 -n 1 --exclusive $*
elif [ "$2" -eq 10000 ]; then
  srun -N 1 -n 16 --exclusive $*
elif [ "$2" -eq 100000 ]; then
  srun -N 1 -n 20 --exclusive $*
elif [ "$2" -eq 1000000 ]; then
  srun -N 2 -n 56 --cpu_bind=sockets --exclusive $*
elif [ "$2" -eq 10000000 ]; then
  srun -N 2 -n 56 --exclusive $*
elif [ "$2" -eq 100000000 ]; then
  srun -N 2 -n 56 --exclusive $*
fi