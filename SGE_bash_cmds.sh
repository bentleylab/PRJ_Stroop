# Julie's potential qalter for all jobs
for i in `qstat -u hoycw | grep ^3 | awk '{print $1}'`; do qalter -l mem_free=3G $i; done

# Check for mem_free status on SGE jobs
for i in `qstat -u hoycw | grep ^3  | awk '{print $1}'`; do echo $i; qstat -j $i  | grep mem_free; done | less
