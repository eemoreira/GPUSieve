#include <bits/stdc++.h>
#define endl '\n'

using namespace std;

const int SEGMENT_SIZE = 1 << 15;

uint64_t segmentedSieve(uint64_t N) {

  const int SQRT = sqrt(N);
  uint64_t s = 2;
  vector<uint8_t> is_prime(SQRT + 1, 1);
  vector<uint8_t> sieve(SEGMENT_SIZE, 1);
  vector<uint32_t> primes;
  primes.reserve(SQRT / log(SQRT));
  is_prime[0] = is_prime[1] = 0;
  uint64_t ans = 0;
  for (uint64_t low = 1; low < N; low += SEGMENT_SIZE) {
    uint64_t high = min(low + SEGMENT_SIZE, N);

    while (s * s < high) {
      if (is_prime[s]) {
        primes.emplace_back(s);
        for (uint64_t j = s * s; j <= SQRT; j += s) {
          is_prime[j] = 0;
        }
      }
      s += 1;
    }

    for (int i = 0; i < SEGMENT_SIZE; i++) sieve[i] = 1;

    for (uint64_t p : primes) {
      uint64_t first = max(2 * p, (low + p - 1) / p * p);
      for (uint64_t i = first; i < high; i += p) {
        sieve[i - low] = 0;
      }
    }

    for (uint64_t i = low; i < high; i++) {
      ans += sieve[i - low];
    }
  }
  return ans - 1;
}

signed main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " n" << endl;
    return 1;
  }

  const uint64_t N = std::atoll(argv[1]);
  cout << segmentedSieve(N) << endl;

  return 0;
}
