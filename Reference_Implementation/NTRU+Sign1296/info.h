#ifndef INFO_H
#define INFO_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include "params.h"

int64_t __rdtsc(void);

#ifdef PARAM_1
// discrete Gaussian
#define K       128
#define KMASK   0x7f
#define w1		13
static const __int128 CDT[w1 + 1] = {
	((__int128)0x1bb040 << 64) ^ (__int128)0x3c4dc2ea906ab1e7,
	((__int128)0x317ba6 << 64) ^ (__int128)0xe28d28c11d21c0f2,
	((__int128)0x3c1cae << 64) ^ (__int128)0x1c63aafdbb01b02b,
	((__int128)0x3f52e0 << 64) ^ (__int128)0x5b191164215da61c,
	((__int128)0x3feccd << 64) ^ (__int128)0xae1320d7dd2d7225,
	((__int128)0x3ffea8 << 64) ^ (__int128)0x78762d0ffa9f467d,
	((__int128)0x3ffff0 << 64) ^ (__int128)0xf66640c0f1f7167d,
	((__int128)0x3fffff << 64) ^ (__int128)0x96e8fbe6fcec75ca,
	((__int128)0x3fffff << 64) ^ (__int128)0xfe37217425bd1884,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffb2e53c2196f17,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffff7ecd89db0b1,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffffff79be6d638,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffffffffa9a4cc0,
	((__int128)0x3fffff << 64) ^ (__int128)0xffffffffffffffff,
};
#define DegG 11
static const uint64_t Coef_G[DegG + 1] = {
	0x00283395212f1fe7,
	0x002ba513b0bf0cdf,
	0x003db2253397e63f,
	0x0004959595cf196f,
	0x00265a8e1133350f,
	0x000461c4deb77207,
	0x0001b75edc624843,
	0x00047b446b5a4491,
	0x0004ae478cf6407f,
	0x003aaab8192ae00f,
	0x003d46835cba7ad7,
	0x1000000000000000 // 2^{60}
};
static const int32_t expDiff_G[DegG] = { 20, 19, 23, 16, 22, 20, 17, 18, 14, 17, 10 };

// rejection sampling
#define DegR 12
static const uint64_t Coef_R[DegR + 1] = {
	0x0f01f2acfe9cd323,
	0x0fe45744e6bb8567,
	0x02fa5c0bf8f06531,
	0x0f9f2a65eb40bd27,
	0x04972cc1c3e51845,
	0x04cb7394652d1f23,
	0x0461c6bdc294a32d,
	0x0dbaf728849e074d,
	0x047b446c244a8911,
	0x04ae478cfb8c16a7,
	0x03aaab8192bf6a9d,
	0x01ea341ae5d3e9d3,
	0x0fffffffffffffff
};
static const int32_t expDiff_R[DegR] = { 41, 43, 38, 42, 40, 40, 38, 41, 39, 39, 39, 34 };

#define vartheta2 60
#define BitGamma 7				// gamma = ceiling[beta2 / (sigma * alpha * ln2)] = 127 < 2^{7}
#define ThetaChat 21
#define ValChat 0x172abd9595	// chat = ValChat / 2^{ThetaChat}
#define ThetaEps 54
#define EpsVal 51857			// epsilon = EpsVal / 2^{ThetaEps}
#endif

#ifdef PARAM_2
// discrete Gaussian
#define K       128
#define KMASK   0x7f
#define w1		14
// 아래 CDT랑 Gaussian용 exp 마무리
static const __int128 CDT[w1 + 1] = {
	((__int128)0x1a08c2 << 64) ^ (__int128)0x604c3fe27ee893da,
	((__int128)0x2f3f46 << 64) ^ (__int128)0x06b1e34ff049124b,
	((__int128)0x3ab8fe << 64) ^ (__int128)0x654fda21b440869f,
	((__int128)0x3ed817 << 64) ^ (__int128)0x8877bf4e21604404,
	((__int128)0x3fd3ad << 64) ^ (__int128)0xd2675e97363d6a3a,
	((__int128)0x3ffb81 << 64) ^ (__int128)0xfba1374915bf7fde,
	((__int128)0x3fffb1 << 64) ^ (__int128)0xa541d1dba472035b,
	((__int128)0x3ffffc << 64) ^ (__int128)0x6d3ec13337a3194e,
	((__int128)0x3fffff << 64) ^ (__int128)0xe427d398013a538c,
	((__int128)0x3fffff << 64) ^ (__int128)0xff6f701a097588d9,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffe0c8c0e024136,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffffb84985259f3,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffffff9282cbe33,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffffffff913b431,
	((__int128)0x3fffff << 64) ^ (__int128)0xffffffffffffffff
};
#define DegG 11
static const uint64_t Coef_G[DegG + 1] = {
	0x00039de6de607651,
	0x00024b6a774a383b,
	0x0001e53d8b8ed795,
	0x001511e36d2e994f,
	0x0019c0ad6a565127,
	0x0000dc13a86bfe81,
	0x000325f3470bd8ff,
	0x000266e436e5a273,
	0x0002ee99958455a1,
	0x00015798ee230265,
	0x000346dc5d638845,
	0x1000000000000000 // 2^{60}
};
static const int32_t expDiff_G[DegG] = { 21, 20, 16, 19, 24, 17, 19, 18, 19, 16, 6 };

// rejection sampling
#define DegR 12
static const uint64_t Coef_R[DegR + 1] = {
	0x024f7f697f05484f,
	0x02dc0a849923ed87,
	0x005029480401f763,
	0x07add9d45d019a15,
	0x05464cbe608c776b,
	0x067058c6b8791005,
	0x01b8280cd17c59f1,
	0x00325f356f827305,
	0x0266e43751e3ac37,
	0x05dd332b0f4ca0a1,
	0x02af31dc461185b3,
	0x01a36e2eb1c432c9,
	0x07ffffffffffffff
};
static const int32_t expDiff_R[DegR] = { 40, 43, 35, 40, 39, 41, 42, 35, 37, 39, 38, 34 };

#define vartheta2 59
#define BitGamma 7 // gamma = ceiling[beta2 / (sigma * alpha * ln2)] = 101 < 2^{7}
#define ThetaChat 20
#define ValChat 0x0000000d89bc6421 // chat = ValChat / 2^{ThetaChat} (Note that ValChat < 2^{36})
#define ThetaEps 52
#define EpsVal 41553 // epsilon = EpsVal / 2^{ThetaEps} (Note that EpsVal < 2^{28})
#endif

#ifdef PARAM_3
// discrete Gaussian
#define K       256
#define KMASK   0xff
#define w1		9
// 아래 CDT랑 Gaussian용 exp 마무리
static const __int128 CDT[w1 + 1] = {
	((__int128)0x22a9cd << 64) ^ (__int128)0x8ba7d055efa8afbb,
	((__int128)0x39234c << 64) ^ (__int128)0x61cf16641dbdf532,
	((__int128)0x3f437c << 64) ^ (__int128)0x2c2e99fb47bc8d28,
	((__int128)0x3ff72b << 64) ^ (__int128)0x199548735b2e2256,
	((__int128)0x3fffd2 << 64) ^ (__int128)0xc0a85a3fb5e4fed4,
	((__int128)0x3fffff << 64) ^ (__int128)0x9de292c9852cb7b6,
	((__int128)0x3fffff << 64) ^ (__int128)0xffa64e7aafbed067,
	((__int128)0x3fffff << 64) ^ (__int128)0xffffdd7d4c74337d,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffffffa6aac9c79,
	((__int128)0x3fffff << 64) ^ (__int128)0xffffffffffffffff,
};
#define DegG 12
static const uint64_t Coef_G[DegG + 1] = {
	0x003e1f87a686ae07,
	0x00289ab9294f64cf,
	0x00022e02f43e5abf,
	0x00065ce0ade4347f,
	0x002115479cbf0a3f,
	0x000131710852e9a9,
	0x00268d07e68e38af,
	0x00215d3fb0774dd7,
	0x0030202c969c9ed7,
	0x003788d0943018df,
	0x003010139a3dbcef,
	0x0037763d2ec76e57,
	0x4000000000000000 // 2^{62}
};
static const int32_t expDiff_G[DegG] = { 22, 25, 19, 18, 25, 15, 20, 19, 19, 19, 18, 9 };

// rejection sampling
#define DegR 13
static const uint64_t Coef_R[DegR + 1] = {
	0x26b957f42ba7de7d,
	0x018883bef8e518c7,
	0x2c53d8ca62a7911b,
	0x11a9b6d418eadb85,
	0x197c855e77b87363,
	0x08458ed2babb5639,
	0x1317195872d85363,
	0x268d08c5ab3b5025,
	0x215d3fb7c9508ad1,
	0x30202c96ea4845f7,
	0x0de234250c7ceb5f,
	0x0c0404e68f709ba5,
	0x037763d2ec76e62b,
	0x3fffffffffffffff
};
static const int32_t expDiff_R[DegR] = { 46, 36, 42, 40, 42, 39, 39, 40, 39, 41, 39, 40, 33 };

#define vartheta2 62
#define BitGamma 8 // gamma = ceiling[beta2 / (sigma * alpha * ln2)] = 164 < 2^{8}
#define ThetaChat 20
#define ValChat 0x1998682d4d // chat = ValChat / 2^{ThetaChat}
#define ThetaEps 56
#define EpsVal 33025 // epsilon = EpsVal / 2^{ThetaEps}

#endif

#endif