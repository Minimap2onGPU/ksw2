#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/sysinfo.h>
#include <omp.h>
#include <hip/hip_runtime.h>

#define ELAPSED(t1,t2) (t2.tv_sec-t1.tv_sec + (t2.tv_usec-t1.tv_usec)*1E-6)

#define HIPCHECK(cmd) \
{\
    hipError_t error  = cmd;\
    if (error != hipSuccess) { \
        fprintf(stderr, "error: '%s'(%d) at %s:%d\n", hipGetErrorString(error), error,__FILE__, __LINE__); \
        exit(EXIT_FAILURE);\
      }\
}

int main (int argc, char *argv[])
{
    int ret, nthread=1;
    int nbuf=1, bufsize_mb=16, niter=100, pinned=1, async_copy=0;
    float avail_mem;
    float *buf, *smbuf;
    float **d_buf;
    size_t bufsize_bytes, nelem;
    struct timeval t1, t2;
    struct sysinfo s;
    hipDeviceProp_t prop;
    hipStream_t stream;

    // Get number of threads
#pragma omp parallel 
    {
        nthread = omp_get_num_threads();
    }
    // ---------------------------------------------------------------
    // Parse arguments
    // ---------------------------------------------------------------
    if (argc < 5) {
        printf ("Usage: %s <nbuf> <bufsize_mb> <pinned> <niter> <async_copy>\n"
                "where,\n"
                "    nbuf       = number of buffers to copy from HtoD (Default: 1)\n"
                "    bufsize_mb = size of each buffer in MB (Default: 16)\n"
                "    pinned     = 0 for pageable host buffer, 1 forpinned host buffer (Default: 1)\n"
                "    niter      = number of iterations of timing loop (Default: 100)\n"
                "    async_copy = 0 for hipMemcpy for transfers, 1 for hipMemcpyAsync (Default: 0)\n",
                argv[0]);
        return 0;
    } 
    nbuf = atoi (argv[1]);
    bufsize_mb = atoi (argv[2]);
    pinned = atoi (argv[3]);
    niter = atoi (argv[4]);
    async_copy = atoi (argv[5]);

    // ---------------------------------------------------------------
    // Error-check arguments 
    // ---------------------------------------------------------------
    if (nbuf<1) {
        printf ("Expecting at least 1 buffer\n");
        return -1;
    }

    // Check buffer size against limits on host and GPU
    // Need to fit 1 buffer in host memory
    ret = sysinfo (&s);
    avail_mem = (float)s.freeram*(float)s.mem_unit/(1048576.f);
    if (bufsize_mb > avail_mem) {
        printf ("Buffer size is too high. Available memory on host = %.1f MB," 
                "requesting %d MB\n", avail_mem, bufsize_mb);
        return -1;
    }
    // Need to fit nbuf buffers in GPU memory
    HIPCHECK (hipGetDeviceProperties(&prop, 0));
    size_t gpumem_avail = prop.totalGlobalMem;
    size_t gpumem_needed = (size_t)nbuf * (size_t)bufsize_mb * (size_t)1048576;
    if (gpumem_avail < gpumem_needed) {
        printf ("gpumem=%zu, needed=%zu\n", gpumem_avail, gpumem_needed);
        printf ("Cannot fit %d buffers of size %d MB in GPU memory\n", 
                nbuf, bufsize_mb);
        return -1;
    }

    if (pinned < 0 || pinned > 1) {
        printf ("pinned must be either 0 or 1\n");
        return -1;
    }
    
    if (niter<1) {
        printf ("Expecting at least 1 iteration\n");
        return -1;
    }

    if (async_copy < 0 || async_copy > 1) {
        printf ("async_copy must be either 0 or 1\n");
        return -1;
    }

    // ---------------------------------------------------------------
    // Allocate a host buffer to copy from, nbuf buffers on the device
    // ---------------------------------------------------------------
    bufsize_bytes = (size_t) bufsize_mb * 1048576;
    nelem = bufsize_bytes/sizeof(float);
    if (async_copy) {
        HIPCHECK (hipStreamCreate (&stream));
    }

    // Allocate and initialize host buffer
    if (!pinned) {
        buf = (float *) malloc (bufsize_bytes);
    } else {
        HIPCHECK (hipHostMalloc (&buf, bufsize_bytes, hipHostMallocDefault));
    }

    // Initialize host buffer
    for (int i=0; i<nelem; i++) {
        buf[i] = (float)i;
    }

    // Allocate device buffer(s)
    d_buf = (float **) malloc (nbuf * sizeof (float *));
    for (int i=0; i<nbuf; i++) {
        HIPCHECK (hipMalloc (&d_buf[i], bufsize_bytes));
    }

    // ---------------------------------------------------------------
    // Time memcpy of multiple small  buffers
    // ---------------------------------------------------------------
    // Amount of data transferred 
    float gbs = ((float)niter*nbuf*(float)bufsize_bytes)/(1024.f*1024.f*1024.f);

    gettimeofday (&t1, NULL);
    for (int iter=0; iter<niter; iter++) {
        #pragma omp parallel for shared(d_buf, buf, bufsize_bytes) num_threads(nthread)
        for (int i=0; i<nbuf; i++) {
            if (!async_copy) {
                HIPCHECK (hipMemcpy (d_buf[i], buf, bufsize_bytes, hipMemcpyHostToDevice));
            } else {
                HIPCHECK (hipMemcpyAsync (d_buf[i], buf, bufsize_bytes, hipMemcpyHostToDevice, stream));
            }
        }
    }
    if (async_copy) {
        hipStreamSynchronize (stream);
    }
    gettimeofday (&t2, NULL);
    printf ("nthread=%d pinned=%d niter=%d async_copy=%d nbuf=%d bufsize=%d MB, bw=%.3f GB/s\n", 
             nthread, pinned, niter, async_copy, nbuf, bufsize_mb, gbs/ELAPSED(t1,t2));

    // ---------------------------------------------------------------
    // Cleanup
    // ---------------------------------------------------------------
    for (int i=0; i<nbuf; i++) {
        HIPCHECK (hipFree (d_buf[i]));
    }
    free (d_buf);

    if (!pinned) {
        free(buf);
    } else {
        HIPCHECK (hipHostFree (buf));
    }
    if (async_copy) {
        HIPCHECK (hipStreamDestroy (stream));
    }

    return 0;
}
