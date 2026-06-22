#include <assert.h>
#include <cuda.h>
#include <math.h>
#include <stdio.h>
#include <sys/select.h>
#include <sys/time.h>
#include <unistd.h>

using ull = unsigned long long;

#define NUM_THREADS_PER_BLOCK 256
#define SQRT_MAX (1 << 20)
const int SEGMENT_SIZE_BYTE = 1 << 10;
const int SEGMENT_SIZE_BIT = 1 << 10;

__device__ int smallPrimes[SQRT_MAX];

__global__ void sieveRangeKernel(ull *sieve_cnt, ull limit, int num_primes) {

  __shared__ uint8_t is_prime[SEGMENT_SIZE_BYTE];

  for (int i = threadIdx.x; i < SEGMENT_SIZE_BYTE; i += blockDim.x) {
    is_prime[i] = 1;
  }

  __syncthreads();

  const ull L = 1ULL * blockIdx.x * SEGMENT_SIZE_BYTE;
  const ull R = min(limit, L + SEGMENT_SIZE_BYTE);
  for (int k = threadIdx.x; k < num_primes; k += blockDim.x) {
    int prime = smallPrimes[k];

    const ull first =
        max((ull)prime * prime, (L + prime - 1) / prime * prime);
    for (ull j = first; j < R; j += prime) {
      is_prime[j - L] = 0;
    }
  }

  __syncthreads();
  ull local_cnt = 0;
  for (int i = threadIdx.x; i < SEGMENT_SIZE_BYTE && L + i < R; i += blockDim.x) {
	  /*
		 print primes
	  if (is_prime[i]) {
		  printf("%llu\n", i + L);
	  }
	  */
      local_cnt += is_prime[i];
  }
  atomicAdd(sieve_cnt, local_cnt);
}

__global__ void sieveRangeBitKernel(ull *sieve_cnt,
                                    ull limit, int num_primes) {

  __shared__ uint32_t is_prime[SEGMENT_SIZE_BIT];

  for (int i = threadIdx.x; i < SEGMENT_SIZE_BIT; i += blockDim.x) {
    is_prime[i] = 0xFFFFFFFF;
  }

  __syncthreads();

  const ull L = 1ULL * blockIdx.x * SEGMENT_SIZE_BIT * 32;
  const ull R = min(limit, L + SEGMENT_SIZE_BIT * 32);
  for (int k = threadIdx.x; k < num_primes; k += blockDim.x) {
    int prime = smallPrimes[k];

    const ull first =
        max((ull)prime * prime, (L + prime - 1) / prime * prime);
    for (ull j = first; j < R; j += prime) {
      int offset = (j - L) >> 5;
      int bit = (j - L) & 31;
      atomicAnd(&is_prime[offset], ~(1U << bit));
    }
  }

  __syncthreads();
  ull local_cnt = 0;
  for (int i = threadIdx.x; i < SEGMENT_SIZE_BIT; i += blockDim.x) {
    ull x = L + (ull)i * 32;
    if (x >= R) {
      break;
    }
    int val = is_prime[i];
    if (x + 31 >= R) {
      int valid_bits = R - x;
      val &= (1U << valid_bits) - 1;
    }
	/*
	   print primes
	for(int b=0;b<32;b++) {
		int cur = (val>>b);
		if(cur&1) printf("%llu\n",x+b);
	}
	*/
    local_cnt += __popc(val);
  }
  atomicAdd(sieve_cnt, local_cnt);
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <limit> <file>\n", argv[0]);
    return 1;
  }

  const ull N = atoll(argv[1]);
  const char *filename = argv[2];
  // compute small primes up to sqrt(N) on the host
  int sqrtN = (int)std::sqrt(N);
  assert(sqrtN < SQRT_MAX);
  bool is_prime[sqrtN + 1];
  for (int i = 0; i <= sqrtN; ++i) {
    is_prime[i] = 1;
  }
  is_prime[0] = is_prime[1] = 0;
  for (int i = 2; i * i <= sqrtN; ++i) {
    if (is_prime[i]) {
      for (int j = i * i; j <= sqrtN; j += i) {
        is_prime[j] = 0;
      }
    }
  }

  int primes[SQRT_MAX];
  memset(primes, 0, sizeof(primes));

  auto byte_sieve = [&](ull limit) {
    timeval start, end;
    gettimeofday(&start, NULL);

    int p = 0;
    for (ull i = 2; i * i < limit; i++) {
      if (is_prime[i]) {
        primes[p++] = i;
      }
    }

    ull *d_sieve_cnt;
    ull sieve_cnt = 0;

    cudaMalloc(&d_sieve_cnt, sizeof(ull));
    cudaMemcpy(d_sieve_cnt, &sieve_cnt, sizeof(ull),
               cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(smallPrimes, primes, p * sizeof(int));
    dim3 dimBlock(NUM_THREADS_PER_BLOCK);
    dim3 dimGrid((limit + SEGMENT_SIZE_BYTE - 1) / SEGMENT_SIZE_BYTE);

    sieveRangeKernel<<<dimGrid, dimBlock>>>(d_sieve_cnt, limit, p);

    cudaError_t err = cudaGetLastError();
    printf("launch: %s\n", cudaGetErrorString(err));

    err = cudaDeviceSynchronize();
    printf("sync: %s\n", cudaGetErrorString(err));

    cudaMemcpy(&sieve_cnt, d_sieve_cnt, sizeof(int), cudaMemcpyDeviceToHost);
    ull ans = (ull)sieve_cnt - 2; // exclude 0 and 1
    printf("Number of primes up to %llu: %llu\n", limit, ans);
    cudaFree(d_sieve_cnt);
    gettimeofday(&end, NULL);
    double elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0 +
                         (end.tv_usec - start.tv_usec) / 1000.0;

    printf("Elapsed time: %.2f ms\n", elapsedTime);
    return elapsedTime;
  };

  auto bit_sieve = [&](ull limit) {
    timeval start, end;
    gettimeofday(&start, NULL);

    ull *d_sieve_cnt;
    ull sieve_cnt = 0;

    int p = 0;
    for (ull i = 2; i * i < limit; i++) {
      if (is_prime[i]) {
        primes[p++] = i;
      }
    }

    cudaMalloc(&d_sieve_cnt, sizeof(ull));
    cudaMemcpy(d_sieve_cnt, &sieve_cnt, sizeof(ull),
               cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(smallPrimes, primes, p * sizeof(int));
    dim3 dimBlock(NUM_THREADS_PER_BLOCK);
    dim3 dimGrid((limit + (SEGMENT_SIZE_BIT * 32) - 1) / (SEGMENT_SIZE_BIT * 32));

    sieveRangeBitKernel<<<dimGrid, dimBlock>>>(d_sieve_cnt, limit, p);

    cudaError_t err = cudaGetLastError();
    printf("launch: %s\n", cudaGetErrorString(err));

    err = cudaDeviceSynchronize();
    printf("sync: %s\n", cudaGetErrorString(err));

    cudaMemcpy(&sieve_cnt, d_sieve_cnt, sizeof(int), cudaMemcpyDeviceToHost);
    ull ans = (ull)sieve_cnt - 2; // exclude 0 and 1
    printf("Number of primes up to %llu: %llu\n", limit, ans);
    cudaFree(d_sieve_cnt);
    gettimeofday(&end, NULL);
    double elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0 +
                         (end.tv_usec - start.tv_usec) / 1000.0;
    printf("Elapsed time: %.2f ms\n", elapsedTime);
    return elapsedTime;
  };

  FILE *f = fopen(filename, "w");
  fprintf(f, "N,Byte Sieve Time (ms),Bit Sieve Time (ms)\n");

  double byte_time = byte_sieve(N);
  double bit_time = bit_sieve(N);
#ifndef SINGLE
  int exp = 2;
  while (true) {
    const ull limit = 1ULL << exp;
    if (limit > N) {
      break;
    }

    printf("Running byte sieve on limit...\n");
	printf("Limit: %llu\n",limit);
    byte_time = byte_sieve(limit);
    printf("====================================\n ");
    printf("Running bit sieve on limit...\n");
	printf("Limit: %llu\n",limit);
    printf("\n\n\n\n");

    bit_time = bit_sieve(limit);
    fprintf(f, "%llu,%.2f,%.2f\n", limit, byte_time, bit_time);
    exp += 2;

  }
#else
  printf("Testing N = %llu\n", N);
  printf("Running byte sieve...\n");
  byte_time = byte_sieve(N);
  printf("====================================\n");
  printf("Running bit sieve...\n");
  bit_time = bit_sieve(N);
  fprintf(f, "%llu,%.2f,%.2f\n", N, byte_time, bit_time);
  printf("\n\n\n\n");
#endif


  return 0;
}
