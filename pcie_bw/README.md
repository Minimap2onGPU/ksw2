# pcie_bw
This standalone helps evaluate the following:
- Is transferring one large buffer faster than transferring 16 small buffers 
  from host to device?
- Is it better to have host buffers in pinned memory even for hipMemcpy calls?
- Is it less efficient if multiple threads launched those small buffer transfers?

## How to build and run

```
make
OMP_NUM_THREADS=<nthread> ./test <nbuf> <bufsize_mb> <pinned> <iter>
```
where,
- nbuf = number of buffers to copy from HtoD (Default: 1)
- bufsize_mb = size of each buffer in MB (Default: 16)
- pinned =0 use pageable host buffer, =1 use pinned host buffer (Default: 1)
- niter  = number of iterations of timing loop (Default: 100)
- async_copy = 0 for hipMemcpy for transfers, 1 for hipMemcpyAsync (Default: 0)

A script is also included to test a range of buffer sizes. To run the script, 
```
./run.sh
```
One could augment the script to vary any or all input arguments as needed.
