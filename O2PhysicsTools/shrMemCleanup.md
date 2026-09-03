# Clearing Stuck Shared Memory in O2Physics (or generally in the machine)

If a process is ungracefully interrupted (e.g., `Ctrl+C` when running `run_all_wagons.sh`, or if any bug kills it unexpectedly), you might want to clean up some shared memory manually. Some commands that proved useful to me were:

1. **Delete abandoned O2 and FairMQ shared memory segments:**
```bash
rm -rf /dev/shm/fmq* /dev/shm/o2*
```

2. **Flush the system caches:**

Sometimes the whole "stage files from NFS to local SSD" or "read a bunch of AO2Ds in consumer" processes take up the whole shrMem due to "FilePages". This should be handled by the OS with no problem at all, but if O2 complains about it anyways, this could help:
```bash
sync
sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
```

3. **Some diagnostics:**
```bash
df -h /dev/shm
numastat -m
numactl --hardware # Checks per-node free memory
```