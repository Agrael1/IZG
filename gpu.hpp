#pragma once
#include "fwd.hpp"
#include <vector>
#include <array>
#include <unordered_map>
#include <bitset>
#include <thread>
#include <future>
#include <mutex>
#include <functional>

#pragma region XMCOLOR
//-------------------------------------------------------------------------------------
// DirectXMath.h -- SIMD C++ Math library
//
// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
//
// http://go.microsoft.com/fwlink/?LinkID=615560
//-------------------------------------------------------------------------------------

#ifndef __cplusplus
#error DirectX Math requires C++
#endif

#define DIRECTX_MATH_VERSION 314

#if defined(_MSC_VER) && (_MSC_VER < 1910)
#error DirectX Math requires Visual C++ 2017 or later.
#endif

#if defined(_MSC_VER) && !defined(_M_ARM) && !defined(_M_ARM64) && !defined(_M_HYBRID_X86_ARM64) && (!_MANAGED) && (!_M_CEE) && (!defined(_M_IX86_FP) || (_M_IX86_FP > 1)) && !defined(_XM_NO_INTRINSICS_) && !defined(_XM_VECTORCALL_)
#define _XM_VECTORCALL_ 1
#endif

#if _XM_VECTORCALL_
#define XM_CALLCONV __vectorcall
#elif defined(__GNUC__)
#define XM_CALLCONV
#else
#define XM_CALLCONV __fastcall
#endif

#ifndef XM_DEPRECATED
#ifdef __GNUC__
#define XM_DEPRECATED __attribute__ ((deprecated))
#else
#define XM_DEPRECATED __declspec(deprecated("This is deprecated and will be removed in a future version."))
#endif
#endif

#if !defined(_XM_AVX2_INTRINSICS_) && defined(__AVX2__) && !defined(_XM_NO_INTRINSICS_)
#define _XM_AVX2_INTRINSICS_
#endif

#if !defined(_XM_FMA3_INTRINSICS_) && defined(_XM_AVX2_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
#define _XM_FMA3_INTRINSICS_
#endif

#if !defined(_XM_F16C_INTRINSICS_) && defined(_XM_AVX2_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
#define _XM_F16C_INTRINSICS_
#endif

#if !defined(_XM_F16C_INTRINSICS_) && defined(__F16C__) && !defined(_XM_NO_INTRINSICS_)
#define _XM_F16C_INTRINSICS_
#endif

#if defined(_XM_FMA3_INTRINSICS_) && !defined(_XM_AVX_INTRINSICS_)
#define _XM_AVX_INTRINSICS_
#endif

#if defined(_XM_F16C_INTRINSICS_) && !defined(_XM_AVX_INTRINSICS_)
#define _XM_AVX_INTRINSICS_
#endif

#if !defined(_XM_AVX_INTRINSICS_) && defined(__AVX__) && !defined(_XM_NO_INTRINSICS_)
#define _XM_AVX_INTRINSICS_
#endif

#if defined(_XM_AVX_INTRINSICS_) && !defined(_XM_SSE4_INTRINSICS_)
#define _XM_SSE4_INTRINSICS_
#endif

#if defined(_XM_SSE4_INTRINSICS_) && !defined(_XM_SSE3_INTRINSICS_)
#define _XM_SSE3_INTRINSICS_
#endif

#if defined(_XM_SSE3_INTRINSICS_) && !defined(_XM_SSE_INTRINSICS_)
#define _XM_SSE_INTRINSICS_
#endif

#if !defined(_XM_ARM_NEON_INTRINSICS_) && !defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
#if (defined(_M_IX86) || defined(_M_X64) || __i386__ || __x86_64__) && !defined(_M_HYBRID_X86_ARM64)
#define _XM_SSE_INTRINSICS_
#elif defined(_M_ARM) || defined(_M_ARM64) || defined(_M_HYBRID_X86_ARM64) || __arm__ || __aarch64__
#define _XM_ARM_NEON_INTRINSICS_
#elif !defined(_XM_NO_INTRINSICS_)
#error DirectX Math does not support this target
#endif
#endif // !_XM_ARM_NEON_INTRINSICS_ && !_XM_SSE_INTRINSICS_ && !_XM_NO_INTRINSICS_

#if !defined(_XM_NO_XMVECTOR_OVERLOADS_) && (defined(__clang__) || defined(__GNUC__))
#define _XM_NO_XMVECTOR_OVERLOADS_
#endif

#pragma warning(push)
#pragma warning(disable:4514 4820)
// C4514/4820: Off by default noise
#include <math.h>
#include <float.h>
#pragma warning(pop)

#ifndef _XM_NO_INTRINSICS_

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4987)
// C4987: Off by default noise
#include <intrin.h>
#pragma warning(pop)
#endif

#if (defined(__clang__) || defined(__GNUC__)) && (__x86_64__ || __i386__)
#include <cpuid.h>
#endif

#ifdef _XM_SSE_INTRINSICS_
#include <xmmintrin.h>
#include <emmintrin.h>

#ifdef _XM_SSE3_INTRINSICS_
#include <pmmintrin.h>
#endif

#ifdef _XM_SSE4_INTRINSICS_
#include <smmintrin.h>
#endif

#ifdef _XM_AVX_INTRINSICS_
#include <immintrin.h>
#endif

#elif defined(_XM_ARM_NEON_INTRINSICS_)
#if defined(_MSC_VER) && (defined(_M_ARM64) || defined(_M_HYBRID_X86_ARM64))
#include <arm64_neon.h>
#else
#include <arm_neon.h>
#endif
#endif
#endif // !_XM_NO_INTRINSICS_

#include <assert.h>

#pragma warning(push)
#pragma warning(disable : 4005 4668)
// C4005/4668: Old header issue
#include <stdint.h>
#pragma warning(pop)

#ifdef __GNUC__
#define XM_ALIGNED_DATA(x) __attribute__ ((aligned(x)))
#define XM_ALIGNED_STRUCT(x) struct __attribute__ ((aligned(x)))
#else
#define XM_ALIGNED_DATA(x) __declspec(align(x))
#define XM_ALIGNED_STRUCT(x) __declspec(align(x)) struct
#endif

/****************************************************************************
 *
 * Conditional intrinsics
 *
 ****************************************************************************/

#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)

#if defined(_XM_NO_MOVNT_)
#define XM_STREAM_PS( p, a ) _mm_store_ps((p), (a))
#define XM256_STREAM_PS( p, a ) _mm256_store_ps((p), (a))
#define XM_SFENCE()
#else
#define XM_STREAM_PS( p, a ) _mm_stream_ps((p), (a))
#define XM256_STREAM_PS( p, a ) _mm256_stream_ps((p), (a))
#define XM_SFENCE() _mm_sfence()
#endif

#if defined(_XM_FMA3_INTRINSICS_)
#define XM_FMADD_PS( a, b, c ) _mm_fmadd_ps((a), (b), (c))
#define XM_FNMADD_PS( a, b, c ) _mm_fnmadd_ps((a), (b), (c))
#else
#define XM_FMADD_PS( a, b, c ) _mm_add_ps(_mm_mul_ps((a), (b)), (c))
#define XM_FNMADD_PS( a, b, c ) _mm_sub_ps((c), _mm_mul_ps((a), (b)))
#endif

#if defined(_XM_AVX_INTRINSICS_) && defined(_XM_FAVOR_INTEL_)
#define XM_PERMUTE_PS( v, c ) _mm_permute_ps((v), c )
#else
#define XM_PERMUTE_PS( v, c ) _mm_shuffle_ps((v), (v), c )
#endif

#endif // _XM_SSE_INTRINSICS_ && !_XM_NO_INTRINSICS_

#if defined(_XM_ARM_NEON_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)

#if defined(__clang__)
#define XM_PREFETCH( a ) __builtin_prefetch(a)
#elif defined(_MSC_VER)
#define XM_PREFETCH( a ) __prefetch(a)
#else
#define XM_PREFETCH( a )
#endif

#endif // _XM_ARM_NEON_INTRINSICS_ && !_XM_NO_INTRINSICS_

namespace DirectX
{
#pragma warning(push)
#pragma warning(disable:4068 4201 4365 4324 4820)
     // C4068: ignore unknown pragmas
     // C4201: nonstandard extension used : nameless struct/union
     // C4365: Off by default noise
     // C4324/4820: padding warnings

//------------------------------------------------------------------------------
#if defined(_XM_NO_INTRINSICS_)
    struct __vector4
    {
        union
        {
            float       vector4_f32[4];
            uint32_t    vector4_u32[4];
        };
    };
#endif // _XM_NO_INTRINSICS_

    //------------------------------------------------------------------------------
    // Vector intrinsic: Four 32 bit floating point components aligned on a 16 byte
    // boundary and mapped to hardware vector registers
#if defined(_XM_SSE_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    typedef __m128 XMVECTOR;
#elif defined(_XM_ARM_NEON_INTRINSICS_) && !defined(_XM_NO_INTRINSICS_)
    typedef float32x4_t XMVECTOR;
#else
    typedef __vector4 XMVECTOR;
#endif

    // Fix-up for (1st-3rd) XMVECTOR parameters that are pass-in-register for x86, ARM, ARM64, and vector call; by reference otherwise
#if ( defined(_M_IX86) || defined(_M_ARM) || defined(_M_ARM64) || _XM_VECTORCALL_ || __i386__ || __arm__ || __aarch64__ ) && !defined(_XM_NO_INTRINSICS_)
    typedef const XMVECTOR FXMVECTOR;
#else
    typedef const XMVECTOR& FXMVECTOR;
#endif

    // Fix-up for (4th) XMVECTOR parameter to pass in-register for ARM, ARM64, and x64 vector call; by reference otherwise
#if ( defined(_M_ARM) || defined(_M_ARM64) || defined(_M_HYBRID_X86_ARM64) || (_XM_VECTORCALL_ && !defined(_M_IX86) ) || __arm__ || __aarch64__ ) && !defined(_XM_NO_INTRINSICS_)
    typedef const XMVECTOR GXMVECTOR;
#else
    typedef const XMVECTOR& GXMVECTOR;
#endif

    // Fix-up for (5th & 6th) XMVECTOR parameter to pass in-register for ARM64 and vector call; by reference otherwise
#if ( defined(_M_ARM64) || defined(_M_HYBRID_X86_ARM64) || _XM_VECTORCALL_ || __aarch64__ ) && !defined(_XM_NO_INTRINSICS_)
    typedef const XMVECTOR HXMVECTOR;
#else
    typedef const XMVECTOR& HXMVECTOR;
#endif

    // Fix-up for (7th+) XMVECTOR parameters to pass by reference
    typedef const XMVECTOR& CXMVECTOR;

    //------------------------------------------------------------------------------
    // Conversion types for constants
    XM_ALIGNED_STRUCT(16) XMVECTORF32
    {
        union
        {
            float f[4];
            XMVECTOR v;
        };

        inline operator XMVECTOR() const noexcept { return v; }
        inline operator const float* () const noexcept { return f; }
#if !defined(_XM_NO_INTRINSICS_) && defined(_XM_SSE_INTRINSICS_)
        inline operator __m128i() const noexcept { return _mm_castps_si128(v); }
        inline operator __m128d() const noexcept { return _mm_castps_pd(v); }
#endif
    };

    XM_ALIGNED_STRUCT(16) XMVECTORI32
    {
        union
        {
            int32_t i[4];
            XMVECTOR v;
        };

        inline operator XMVECTOR() const noexcept { return v; }
#if !defined(_XM_NO_INTRINSICS_) && defined(_XM_SSE_INTRINSICS_)
        inline operator __m128i() const noexcept { return _mm_castps_si128(v); }
        inline operator __m128d() const noexcept { return _mm_castps_pd(v); }
#endif
    };

    XM_ALIGNED_STRUCT(16) XMVECTORU32
    {
        union
        {
            uint32_t u[4];
            XMVECTOR v;
        };

        inline operator XMVECTOR() const noexcept { return v; }
#if !defined(_XM_NO_INTRINSICS_) && defined(_XM_SSE_INTRINSICS_)
        inline operator __m128i() const noexcept { return _mm_castps_si128(v); }
        inline operator __m128d() const noexcept { return _mm_castps_pd(v); }
#endif
    };

    //------------------------------------------------------------------------------
    // Matrix type: Sixteen 32 bit floating point components aligned on a
    // 16 byte boundary and mapped to four hardware vector registers

    struct XMMATRIX;

    // Fix-up for (1st) XMMATRIX parameter to pass in-register for ARM64 and vector call; by reference otherwise
#if ( defined(_M_ARM64) || defined(_M_HYBRID_X86_ARM64) || _XM_VECTORCALL_ || __aarch64__ ) && !defined(_XM_NO_INTRINSICS_)
    typedef const XMMATRIX FXMMATRIX;
#else
    typedef const XMMATRIX& FXMMATRIX;
#endif

    // Fix-up for (2nd+) XMMATRIX parameters to pass by reference
    typedef const XMMATRIX& CXMMATRIX;

#ifdef _XM_NO_INTRINSICS_
    struct XMMATRIX
#else
    XM_ALIGNED_STRUCT(16) XMMATRIX
#endif
    {
#ifdef _XM_NO_INTRINSICS_
        union
        {
            XMVECTOR r[4];
            struct
            {
                float _11, _12, _13, _14;
                float _21, _22, _23, _24;
                float _31, _32, _33, _34;
                float _41, _42, _43, _44;
            };
            float m[4][4];
        };
#else
        XMVECTOR r[4];
#endif

        XMMATRIX() = default;

        XMMATRIX(const XMMATRIX&) = default;

#if defined(_MSC_VER) && (_MSC_FULL_VER < 191426431)
        XMMATRIX& operator= (const XMMATRIX& M) noexcept { r[0] = M.r[0]; r[1] = M.r[1]; r[2] = M.r[2]; r[3] = M.r[3]; return *this; }
#else
        XMMATRIX& operator=(const XMMATRIX&) = default;

        XMMATRIX(XMMATRIX&&) = default;
        XMMATRIX& operator=(XMMATRIX&&) = default;
#endif
        constexpr XMMATRIX(FXMVECTOR R0, FXMVECTOR R1, FXMVECTOR R2, CXMVECTOR R3) noexcept : r{ R0,R1,R2,R3 } {}
#ifdef _XM_NO_INTRINSICS_
        float       operator() (size_t Row, size_t Column) const noexcept { return m[Row][Column]; }
        float& operator() (size_t Row, size_t Column) noexcept { return m[Row][Column]; }
#endif
        XMMATRIX    XM_CALLCONV     operator* (FXMMATRIX M) const noexcept;
    };

    //------------------------------------------------------------------------------
    // 4D Vector; 32 bit floating point components
    struct XMFLOAT4
    {
        float x;
        float y;
        float z;
        float w;

        XMFLOAT4() = default;

        XMFLOAT4(const XMFLOAT4&) = default;
        XMFLOAT4& operator=(const XMFLOAT4&) = default;

        XMFLOAT4(XMFLOAT4&&) = default;
        XMFLOAT4& operator=(XMFLOAT4&&) = default;

        constexpr XMFLOAT4(float _x, float _y, float _z, float _w) noexcept : x(_x), y(_y), z(_z), w(_w) {}
        explicit XMFLOAT4(const float* pArray) noexcept : x(pArray[0]), y(pArray[1]), z(pArray[2]), w(pArray[3]) {}
    };

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-anonymous-struct"
#pragma clang diagnostic ignored "-Wnested-anon-types"
#endif

    //------------------------------------------------------------------------------
    // 4x4 Matrix: 32 bit floating point components
    struct XMFLOAT4X4
    {
        union
        {
            struct
            {
                float _11, _12, _13, _14;
                float _21, _22, _23, _24;
                float _31, _32, _33, _34;
                float _41, _42, _43, _44;
            };
            float m[4][4];
        };

        XMFLOAT4X4() = default;

        XMFLOAT4X4(const XMFLOAT4X4&) = default;
        XMFLOAT4X4& operator=(const XMFLOAT4X4&) = default;

        XMFLOAT4X4(XMFLOAT4X4&&) = default;
        XMFLOAT4X4& operator=(XMFLOAT4X4&&) = default;

        constexpr XMFLOAT4X4(float m00, float m01, float m02, float m03,
            float m10, float m11, float m12, float m13,
            float m20, float m21, float m22, float m23,
            float m30, float m31, float m32, float m33) noexcept
            : _11(m00), _12(m01), _13(m02), _14(m03),
            _21(m10), _22(m11), _23(m12), _24(m13),
            _31(m20), _32(m21), _33(m22), _34(m23),
            _41(m30), _42(m31), _43(m32), _44(m33) {}
        explicit XMFLOAT4X4(const float* pArray) noexcept;

        float       operator() (size_t Row, size_t Column) const noexcept { return m[Row][Column]; }
        float& operator() (size_t Row, size_t Column) noexcept { return m[Row][Column]; }
    };

    ////////////////////////////////////////////////////////////////////////////////

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#pragma warning(pop)

#ifndef XMGLOBALCONST
#if defined(__GNUC__) && !defined(__MINGW32__)
#define XMGLOBALCONST extern const __attribute__((weak))
#else
#define XMGLOBALCONST extern const __declspec(selectany)
#endif
#endif
    XMGLOBALCONST XMVECTORU32 g_XMMask3 = { { { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000 } } };
    XMGLOBALCONST XMVECTORF32 g_XMOne = { { { 1.0f, 1.0f, 1.0f, 1.0f } } };
    XMGLOBALCONST XMVECTORF32 g_XMOne3 = { { { 1.0f, 1.0f, 1.0f, 0.0f } } };
    XMGLOBALCONST XMVECTORF32 g_XMZero = { { { 0.0f, 0.0f, 0.0f, 0.0f } } };
    
    XMGLOBALCONST XMVECTORI32 g_XMInfinity = { { { 0x7F800000, 0x7F800000, 0x7F800000, 0x7F800000 } } };
    XMGLOBALCONST XMVECTORI32 g_XMQNaN = { { { 0x7FC00000, 0x7FC00000, 0x7FC00000, 0x7FC00000 } } };

    XMGLOBALCONST XMVECTORF32 g_UByteMax = { { { 255.0f, 255.0f, 255.0f, 255.0f } } };

#pragma warning(push)
#pragma warning(disable:4068 4214 4204 4365 4616 4640 6001 6101)
     // C4068/4616: ignore unknown pragmas
     // C4214/4204: nonstandard extension used
     // C4365/4640: Off by default noise
     // C6001/6101: False positives

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wundefined-reinterpret-cast"
#endif
    inline void XM_CALLCONV XMStoreFloat(float* pDestination, FXMVECTOR V) noexcept
    {
        assert(pDestination);
#if defined(_XM_NO_INTRINSICS_)
        *pDestination = XMVectorGetX(V);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        vst1q_lane_f32(pDestination, V, 0);
#elif defined(_XM_SSE_INTRINSICS_)
        _mm_store_ss(pDestination, V);
#endif
    }
    inline void XM_CALLCONV XMStoreFloat4(XMFLOAT4* pDestination, FXMVECTOR V) noexcept
    {
        assert(pDestination);
#if defined(_XM_NO_INTRINSICS_)
        pDestination->x = V.vector4_f32[0];
        pDestination->y = V.vector4_f32[1];
        pDestination->z = V.vector4_f32[2];
        pDestination->w = V.vector4_f32[3];
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        vst1q_f32(reinterpret_cast<float*>(pDestination), V);
#elif defined(_XM_SSE_INTRINSICS_)
        _mm_storeu_ps(&pDestination->x, V);
#endif
    }
    inline XMMATRIX XM_CALLCONV XMLoadFloat4x4(const XMFLOAT4X4* pSource) noexcept
    {
        assert(pSource);
#if defined(_XM_NO_INTRINSICS_)

        XMMATRIX M;
        M.r[0].vector4_f32[0] = pSource->m[0][0];
        M.r[0].vector4_f32[1] = pSource->m[0][1];
        M.r[0].vector4_f32[2] = pSource->m[0][2];
        M.r[0].vector4_f32[3] = pSource->m[0][3];

        M.r[1].vector4_f32[0] = pSource->m[1][0];
        M.r[1].vector4_f32[1] = pSource->m[1][1];
        M.r[1].vector4_f32[2] = pSource->m[1][2];
        M.r[1].vector4_f32[3] = pSource->m[1][3];

        M.r[2].vector4_f32[0] = pSource->m[2][0];
        M.r[2].vector4_f32[1] = pSource->m[2][1];
        M.r[2].vector4_f32[2] = pSource->m[2][2];
        M.r[2].vector4_f32[3] = pSource->m[2][3];

        M.r[3].vector4_f32[0] = pSource->m[3][0];
        M.r[3].vector4_f32[1] = pSource->m[3][1];
        M.r[3].vector4_f32[2] = pSource->m[3][2];
        M.r[3].vector4_f32[3] = pSource->m[3][3];
        return M;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        XMMATRIX M;
        M.r[0] = vld1q_f32(reinterpret_cast<const float*>(&pSource->_11));
        M.r[1] = vld1q_f32(reinterpret_cast<const float*>(&pSource->_21));
        M.r[2] = vld1q_f32(reinterpret_cast<const float*>(&pSource->_31));
        M.r[3] = vld1q_f32(reinterpret_cast<const float*>(&pSource->_41));
        return M;
#elif defined(_XM_SSE_INTRINSICS_)
        XMMATRIX M;
        M.r[0] = _mm_loadu_ps(&pSource->_11);
        M.r[1] = _mm_loadu_ps(&pSource->_21);
        M.r[2] = _mm_loadu_ps(&pSource->_31);
        M.r[3] = _mm_loadu_ps(&pSource->_41);
        return M;
#endif
    }
    inline XMVECTOR XM_CALLCONV XMLoadFloat4(const XMFLOAT4* pSource) noexcept
    {
        assert(pSource);
#if defined(_XM_NO_INTRINSICS_)
        XMVECTOR V;
        V.vector4_f32[0] = pSource->x;
        V.vector4_f32[1] = pSource->y;
        V.vector4_f32[2] = pSource->z;
        V.vector4_f32[3] = pSource->w;
        return V;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vld1q_f32(reinterpret_cast<const float*>(pSource));
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_loadu_ps(&pSource->x);
#endif
    }
    inline bool XM_CALLCONV XMVector4LessOrEqual(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        return (((V1.vector4_f32[0] <= V2.vector4_f32[0]) && (V1.vector4_f32[1] <= V2.vector4_f32[1]) && (V1.vector4_f32[2] <= V2.vector4_f32[2]) && (V1.vector4_f32[3] <= V2.vector4_f32[3])) != 0);
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        uint32x4_t vResult = vcleq_f32(V1, V2);
        uint8x8x2_t vTemp = vzip_u8(vget_low_u8(vResult), vget_high_u8(vResult));
        uint16x4x2_t vTemp2 = vzip_u16(vTemp.val[0], vTemp.val[1]);
        return (vget_lane_u32(vTemp2.val[1], 1) == 0xFFFFFFFFU);
#elif defined(_XM_SSE_INTRINSICS_)
        XMVECTOR vTemp = _mm_cmple_ps(V1, V2);
        return ((_mm_movemask_ps(vTemp) == 0x0f) != 0);
#else
        return XMComparisonAllTrue(XMVector4GreaterOrEqualR(V2, V1));
#endif
}


    //--------------------Vector------------------------------------------------
    inline XMVECTOR XM_CALLCONV XMVectorZero() noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 vResult = { { { 0.0f, 0.0f, 0.0f, 0.0f } } };
        return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vdupq_n_f32(0);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_setzero_ps();
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorAdd(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORF32 Result = { { {
                V1.vector4_f32[0] + V2.vector4_f32[0],
                V1.vector4_f32[1] + V2.vector4_f32[1],
                V1.vector4_f32[2] + V2.vector4_f32[2],
                V1.vector4_f32[3] + V2.vector4_f32[3]
            } } };
        return Result.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vaddq_f32(V1, V2);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_add_ps(V1, V2);
#endif
    }

    inline XMVECTOR XM_CALLCONV XMVectorMergeZW(FXMVECTOR V1, FXMVECTOR V2)
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORU32 Result = { { {
                V1.vector4_u32[2],
                V2.vector4_u32[2],
                V1.vector4_u32[3],
                V2.vector4_u32[3]
            } } };
        return Result.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vzipq_f32(V1, V2).val[1];
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_unpackhi_ps(V1, V2);
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorSet(float x, float y, float z, float w) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 vResult = { { { x, y, z, w } } };
        return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        float32x2_t V0 = vcreate_f32(
            static_cast<uint64_t>(*reinterpret_cast<const uint32_t*>(&x))
            | (static_cast<uint64_t>(*reinterpret_cast<const uint32_t*>(&y)) << 32));
        float32x2_t V1 = vcreate_f32(
            static_cast<uint64_t>(*reinterpret_cast<const uint32_t*>(&z))
            | (static_cast<uint64_t>(*reinterpret_cast<const uint32_t*>(&w)) << 32));
        return vcombine_f32(V0, V1);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_set_ps(w, z, y, x);
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorReplicate(float Value) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 vResult;
        vResult.f[0] =
            vResult.f[1] =
            vResult.f[2] =
            vResult.f[3] = Value;
        return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vdupq_n_f32(Value);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_set_ps1(Value);
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorSplatW(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 vResult;
        vResult.f[0] =
            vResult.f[1] =
            vResult.f[2] =
            vResult.f[3] = V.vector4_f32[3];
        return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vdupq_lane_f32(vget_high_f32(V), 1);
#elif defined(_XM_SSE_INTRINSICS_)
        return XM_PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3));
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorLess(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORU32 Control = { { {
                (V1.vector4_f32[0] < V2.vector4_f32[0]) ? 0xFFFFFFFF : 0,
                (V1.vector4_f32[1] < V2.vector4_f32[1]) ? 0xFFFFFFFF : 0,
                (V1.vector4_f32[2] < V2.vector4_f32[2]) ? 0xFFFFFFFF : 0,
                (V1.vector4_f32[3] < V2.vector4_f32[3]) ? 0xFFFFFFFF : 0
            } } };
        return Control.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vcltq_f32(V1, V2);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_cmplt_ps(V1, V2);
#endif
}
    inline XMVECTOR XM_CALLCONV XMVectorAndInt(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORU32 Result = { { {
                V1.vector4_u32[0] & V2.vector4_u32[0],
                V1.vector4_u32[1] & V2.vector4_u32[1],
                V1.vector4_u32[2] & V2.vector4_u32[2],
                V1.vector4_u32[3] & V2.vector4_u32[3]
            } } };
        return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vandq_u32(V1, V2);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_and_ps(V1, V2);
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorNegate(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORF32 Result = { { {
                -V.vector4_f32[0],
                -V.vector4_f32[1],
                -V.vector4_f32[2],
                -V.vector4_f32[3]
            } } };
        return Result.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vnegq_f32(V);
#elif defined(_XM_SSE_INTRINSICS_)
        XMVECTOR Z;

        Z = _mm_setzero_ps();

        return _mm_sub_ps(Z, V);
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorReciprocalEst(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 Result = { { {
                1.f / V.vector4_f32[0],
                1.f / V.vector4_f32[1],
                1.f / V.vector4_f32[2],
                1.f / V.vector4_f32[3]
            } } };
        return Result.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vrecpeq_f32(V);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_rcp_ps(V);
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorReciprocal(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 Result = { { {
                1.f / V.vector4_f32[0],
                1.f / V.vector4_f32[1],
                1.f / V.vector4_f32[2],
                1.f / V.vector4_f32[3]
            } } };
        return Result.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
#if defined(_M_ARM64) || defined(_M_HYBRID_X86_ARM64) || __aarch64__
        float32x4_t one = vdupq_n_f32(1.0f);
        return vdivq_f32(one, V);
#else
        // 2 iterations of Newton-Raphson refinement
        float32x4_t Reciprocal = vrecpeq_f32(V);
        float32x4_t S = vrecpsq_f32(Reciprocal, V);
        Reciprocal = vmulq_f32(S, Reciprocal);
        S = vrecpsq_f32(Reciprocal, V);
        return vmulq_f32(S, Reciprocal);
#endif
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_div_ps(g_XMOne, V);
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorSplatX(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 vResult;
        vResult.f[0] =
            vResult.f[1] =
            vResult.f[2] =
            vResult.f[3] = V.vector4_f32[0];
        return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vdupq_lane_f32(vget_low_f32(V), 0);
#elif defined(_XM_AVX2_INTRINSICS_) && defined(_XM_FAVOR_INTEL_)
        return _mm_broadcastss_ps(V);
#elif defined(_XM_SSE_INTRINSICS_)
        return XM_PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0));
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorSplatY(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 vResult;
        vResult.f[0] =
            vResult.f[1] =
            vResult.f[2] =
            vResult.f[3] = V.vector4_f32[1];
        return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vdupq_lane_f32(vget_low_f32(V), 1);
#elif defined(_XM_SSE_INTRINSICS_)
        return XM_PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1));
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorSplatZ(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 vResult;
        vResult.f[0] =
            vResult.f[1] =
            vResult.f[2] =
            vResult.f[3] = V.vector4_f32[2];
        return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vdupq_lane_f32(vget_high_f32(V), 0);
#elif defined(_XM_SSE_INTRINSICS_)
        return XM_PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2));
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorMultiply(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 Result = { { {
                V1.vector4_f32[0] * V2.vector4_f32[0],
                V1.vector4_f32[1] * V2.vector4_f32[1],
                V1.vector4_f32[2] * V2.vector4_f32[2],
                V1.vector4_f32[3] * V2.vector4_f32[3]
            } } };
        return Result.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vmulq_f32(V1, V2);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_mul_ps(V1, V2);
#endif
}
    inline XMVECTOR XM_CALLCONV XMVectorMultiplyAdd(FXMVECTOR V1, FXMVECTOR V2, FXMVECTOR V3) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 Result = { { {
                V1.vector4_f32[0] * V2.vector4_f32[0] + V3.vector4_f32[0],
                V1.vector4_f32[1] * V2.vector4_f32[1] + V3.vector4_f32[1],
                V1.vector4_f32[2] * V2.vector4_f32[2] + V3.vector4_f32[2],
                V1.vector4_f32[3] * V2.vector4_f32[3] + V3.vector4_f32[3]
            } } };
        return Result.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
#if defined(_M_ARM64) || defined(_M_HYBRID_X86_ARM64) || __aarch64__
        return vfmaq_f32(V3, V1, V2);
#else
        return vmlaq_f32(V3, V1, V2);
#endif
#elif defined(_XM_SSE_INTRINSICS_)
        return XM_FMADD_PS(V1, V2, V3);
#endif
}
    inline XMVECTOR XM_CALLCONV XMVectorNegativeMultiplySubtract(FXMVECTOR V1, FXMVECTOR V2, FXMVECTOR V3) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 Result = { { {
                V3.vector4_f32[0] - (V1.vector4_f32[0] * V2.vector4_f32[0]),
                V3.vector4_f32[1] - (V1.vector4_f32[1] * V2.vector4_f32[1]),
                V3.vector4_f32[2] - (V1.vector4_f32[2] * V2.vector4_f32[2]),
                V3.vector4_f32[3] - (V1.vector4_f32[3] * V2.vector4_f32[3])
            } } };
        return Result;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
#if defined(_M_ARM64) || defined(_M_HYBRID_X86_ARM64) || __aarch64__
        return vfmsq_f32(V3, V1, V2);
#else
        return vmlsq_f32(V3, V1, V2);
#endif
#elif defined(_XM_SSE_INTRINSICS_)
        return XM_FNMADD_PS(V1, V2, V3);
#endif
}
    inline XMVECTOR XM_CALLCONV XMVectorSubtract(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORF32 Result = { { {
                V1.vector4_f32[0] - V2.vector4_f32[0],
                V1.vector4_f32[1] - V2.vector4_f32[1],
                V1.vector4_f32[2] - V2.vector4_f32[2],
                V1.vector4_f32[3] - V2.vector4_f32[3]
            } } };
        return Result.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vsubq_f32(V1, V2);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_sub_ps(V1, V2);
#endif
}
    inline XMVECTOR XM_CALLCONV XMVectorLerp(FXMVECTOR V0, FXMVECTOR V1, float t) noexcept
    {
        // V0 + t * (V1 - V0)

#if defined(_XM_NO_INTRINSICS_)

        XMVECTOR Scale = XMVectorReplicate(t);
        XMVECTOR Length = XMVectorSubtract(V1, V0);
        return XMVectorMultiplyAdd(Length, Scale, V0);

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        XMVECTOR L = vsubq_f32(V1, V0);
        return vmlaq_n_f32(V0, L, t);
#elif defined(_XM_SSE_INTRINSICS_)
        XMVECTOR L = _mm_sub_ps(V1, V0);
        XMVECTOR S = _mm_set_ps1(t);
        return XM_FMADD_PS(L, S, V0);
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorEqual(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORU32 Control = { { {
                (V1.vector4_f32[0] == V2.vector4_f32[0]) ? 0xFFFFFFFF : 0,
                (V1.vector4_f32[1] == V2.vector4_f32[1]) ? 0xFFFFFFFF : 0,
                (V1.vector4_f32[2] == V2.vector4_f32[2]) ? 0xFFFFFFFF : 0,
                (V1.vector4_f32[3] == V2.vector4_f32[3]) ? 0xFFFFFFFF : 0,
            } } };
        return Control.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vceqq_f32(V1, V2);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_cmpeq_ps(V1, V2);
#endif
}
    inline XMVECTOR XM_CALLCONV XMVectorEqualInt(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORU32 Control = { { {
                (V1.vector4_u32[0] == V2.vector4_u32[0]) ? 0xFFFFFFFF : 0,
                (V1.vector4_u32[1] == V2.vector4_u32[1]) ? 0xFFFFFFFF : 0,
                (V1.vector4_u32[2] == V2.vector4_u32[2]) ? 0xFFFFFFFF : 0,
                (V1.vector4_u32[3] == V2.vector4_u32[3]) ? 0xFFFFFFFF : 0,
            } } };
        return Control.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vceqq_u32(V1, V2);
#elif defined(_XM_SSE_INTRINSICS_)
        __m128i V = _mm_cmpeq_epi32(_mm_castps_si128(V1), _mm_castps_si128(V2));
        return _mm_castsi128_ps(V);
#endif
}
    inline XMVECTOR XM_CALLCONV XMVectorSelect(FXMVECTOR V1, FXMVECTOR V2, FXMVECTOR Control) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORU32 Result = { { {
                (V1.vector4_u32[0] & ~Control.vector4_u32[0]) | (V2.vector4_u32[0] & Control.vector4_u32[0]),
                (V1.vector4_u32[1] & ~Control.vector4_u32[1]) | (V2.vector4_u32[1] & Control.vector4_u32[1]),
                (V1.vector4_u32[2] & ~Control.vector4_u32[2]) | (V2.vector4_u32[2] & Control.vector4_u32[2]),
                (V1.vector4_u32[3] & ~Control.vector4_u32[3]) | (V2.vector4_u32[3] & Control.vector4_u32[3]),
            } } };
        return Result.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vbslq_f32(Control, V2, V1);
#elif defined(_XM_SSE_INTRINSICS_)
        XMVECTOR vTemp1 = _mm_andnot_ps(Control, V1);
        XMVECTOR vTemp2 = _mm_and_ps(V2, Control);
        return _mm_or_ps(vTemp1, vTemp2);
#endif
}
    inline XMVECTOR XM_CALLCONV XMVectorSqrt(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMVECTORF32 Result = { { {
                sqrtf(V.vector4_f32[0]),
                sqrtf(V.vector4_f32[1]),
                sqrtf(V.vector4_f32[2]),
                sqrtf(V.vector4_f32[3])
            } } };
        return Result.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        // 3 iterations of Newton-Raphson refinment of sqrt
        float32x4_t S0 = vrsqrteq_f32(V);
        float32x4_t P0 = vmulq_f32(V, S0);
        float32x4_t R0 = vrsqrtsq_f32(P0, S0);
        float32x4_t S1 = vmulq_f32(S0, R0);
        float32x4_t P1 = vmulq_f32(V, S1);
        float32x4_t R1 = vrsqrtsq_f32(P1, S1);
        float32x4_t S2 = vmulq_f32(S1, R1);
        float32x4_t P2 = vmulq_f32(V, S2);
        float32x4_t R2 = vrsqrtsq_f32(P2, S2);
        float32x4_t S3 = vmulq_f32(S2, R2);

        XMVECTOR VEqualsInfinity = XMVectorEqualInt(V, g_XMInfinity.v);
        XMVECTOR VEqualsZero = XMVectorEqual(V, vdupq_n_f32(0));
        XMVECTOR Result = vmulq_f32(V, S3);
        XMVECTOR Select = XMVectorEqualInt(VEqualsInfinity, VEqualsZero);
        return XMVectorSelect(V, Result, Select);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_sqrt_ps(V);
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorMax(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORF32 Result = { { {
                (V1.vector4_f32[0] > V2.vector4_f32[0]) ? V1.vector4_f32[0] : V2.vector4_f32[0],
                (V1.vector4_f32[1] > V2.vector4_f32[1]) ? V1.vector4_f32[1] : V2.vector4_f32[1],
                (V1.vector4_f32[2] > V2.vector4_f32[2]) ? V1.vector4_f32[2] : V2.vector4_f32[2],
                (V1.vector4_f32[3] > V2.vector4_f32[3]) ? V1.vector4_f32[3] : V2.vector4_f32[3]
            } } };
        return Result.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        return vmaxq_f32(V1, V2);
#elif defined(_XM_SSE_INTRINSICS_)
        return _mm_max_ps(V1, V2);
#endif
}
    inline XMVECTOR XM_CALLCONV XMVectorPow(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTORF32 Result = { { {
                powf(V1.vector4_f32[0], V2.vector4_f32[0]),
                powf(V1.vector4_f32[1], V2.vector4_f32[1]),
                powf(V1.vector4_f32[2], V2.vector4_f32[2]),
                powf(V1.vector4_f32[3], V2.vector4_f32[3])
            } } };
        return Result.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        XMVECTORF32 vResult = { { {
                powf(vgetq_lane_f32(V1, 0), vgetq_lane_f32(V2, 0)),
                powf(vgetq_lane_f32(V1, 1), vgetq_lane_f32(V2, 1)),
                powf(vgetq_lane_f32(V1, 2), vgetq_lane_f32(V2, 2)),
                powf(vgetq_lane_f32(V1, 3), vgetq_lane_f32(V2, 3))
            } } };
        return vResult.v;
#elif defined(_XM_SSE_INTRINSICS_)
        XM_ALIGNED_DATA(16) float a[4];
        XM_ALIGNED_DATA(16) float b[4];
        _mm_store_ps(a, V1);
        _mm_store_ps(b, V2);
        XMVECTOR vResult = _mm_setr_ps(
            powf(a[0], b[0]),
            powf(a[1], b[1]),
            powf(a[2], b[2]),
            powf(a[3], b[3]));
        return vResult;
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorClamp(FXMVECTOR V, FXMVECTOR Min, FXMVECTOR Max) noexcept
    {
        assert(XMVector4LessOrEqual(Min, Max));

#if defined(_XM_NO_INTRINSICS_)

        XMVECTOR Result;
        Result = XMVectorMax(Min, V);
        Result = XMVectorMin(Max, Result);
        return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        XMVECTOR vResult;
        vResult = vmaxq_f32(Min, V);
        vResult = vminq_f32(Max, vResult);
        return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
        XMVECTOR vResult;
        vResult = _mm_max_ps(Min, V);
        vResult = _mm_min_ps(Max, vResult);
        return vResult;
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVectorSaturate(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        const XMVECTOR Zero = XMVectorZero();

        return XMVectorClamp(V, Zero, g_XMOne.v);

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        // Set <0 to 0
        XMVECTOR vResult = vmaxq_f32(V, vdupq_n_f32(0));
        // Set>1 to 1
        return vminq_f32(vResult, vdupq_n_f32(1.0f));
#elif defined(_XM_SSE_INTRINSICS_)
        // Set <0 to 0
        XMVECTOR vResult = _mm_max_ps(V, g_XMZero);
        // Set>1 to 1
        return _mm_min_ps(vResult, g_XMOne);
#endif
    }

    //--------------------VectorX-----------------------------------------------
    inline XMVECTOR XM_CALLCONV XMVector4Transform(FXMVECTOR V, FXMMATRIX M) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        float fX = (M.m[0][0] * V.vector4_f32[0]) + (M.m[1][0] * V.vector4_f32[1]) + (M.m[2][0] * V.vector4_f32[2]) + (M.m[3][0] * V.vector4_f32[3]);
        float fY = (M.m[0][1] * V.vector4_f32[0]) + (M.m[1][1] * V.vector4_f32[1]) + (M.m[2][1] * V.vector4_f32[2]) + (M.m[3][1] * V.vector4_f32[3]);
        float fZ = (M.m[0][2] * V.vector4_f32[0]) + (M.m[1][2] * V.vector4_f32[1]) + (M.m[2][2] * V.vector4_f32[2]) + (M.m[3][2] * V.vector4_f32[3]);
        float fW = (M.m[0][3] * V.vector4_f32[0]) + (M.m[1][3] * V.vector4_f32[1]) + (M.m[2][3] * V.vector4_f32[2]) + (M.m[3][3] * V.vector4_f32[3]);
        XMVECTORF32 vResult = { { { fX, fY, fZ, fW } } };
        return vResult.v;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        float32x2_t VL = vget_low_f32(V);
        XMVECTOR vResult = vmulq_lane_f32(M.r[0], VL, 0); // X
        vResult = vmlaq_lane_f32(vResult, M.r[1], VL, 1); // Y
        float32x2_t VH = vget_high_f32(V);
        vResult = vmlaq_lane_f32(vResult, M.r[2], VH, 0); // Z
        return vmlaq_lane_f32(vResult, M.r[3], VH, 1); // W
#elif defined(_XM_SSE_INTRINSICS_)
        XMVECTOR vResult = XM_PERMUTE_PS(V, _MM_SHUFFLE(3, 3, 3, 3)); // W
        vResult = _mm_mul_ps(vResult, M.r[3]);
        XMVECTOR vTemp = XM_PERMUTE_PS(V, _MM_SHUFFLE(2, 2, 2, 2)); // Z
        vResult = XM_FMADD_PS(vTemp, M.r[2], vResult);
        vTemp = XM_PERMUTE_PS(V, _MM_SHUFFLE(1, 1, 1, 1)); // Y
        vResult = XM_FMADD_PS(vTemp, M.r[1], vResult);
        vTemp = XM_PERMUTE_PS(V, _MM_SHUFFLE(0, 0, 0, 0)); // X
        vResult = XM_FMADD_PS(vTemp, M.r[0], vResult);
        return vResult;
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVector3Dot(FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        float fValue = V1.vector4_f32[0] * V2.vector4_f32[0] + V1.vector4_f32[1] * V2.vector4_f32[1] + V1.vector4_f32[2] * V2.vector4_f32[2];
        XMVECTORF32 vResult;
        vResult.f[0] =
            vResult.f[1] =
            vResult.f[2] =
            vResult.f[3] = fValue;
        return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        float32x4_t vTemp = vmulq_f32(V1, V2);
        float32x2_t v1 = vget_low_f32(vTemp);
        float32x2_t v2 = vget_high_f32(vTemp);
        v1 = vpadd_f32(v1, v1);
        v2 = vdup_lane_f32(v2, 0);
        v1 = vadd_f32(v1, v2);
        return vcombine_f32(v1, v1);
#elif defined(_XM_SSE4_INTRINSICS_)
        return _mm_dp_ps(V1, V2, 0x7f);
#elif defined(_XM_SSE3_INTRINSICS_)
        XMVECTOR vTemp = _mm_mul_ps(V1, V2);
        vTemp = _mm_and_ps(vTemp, g_XMMask3);
        vTemp = _mm_hadd_ps(vTemp, vTemp);
        return _mm_hadd_ps(vTemp, vTemp);
#elif defined(_XM_SSE_INTRINSICS_)
        // Perform the dot product
        XMVECTOR vDot = _mm_mul_ps(V1, V2);
        // x=Dot.vector4_f32[1], y=Dot.vector4_f32[2]
        XMVECTOR vTemp = XM_PERMUTE_PS(vDot, _MM_SHUFFLE(2, 1, 2, 1));
        // Result.vector4_f32[0] = x+y
        vDot = _mm_add_ss(vDot, vTemp);
        // x=Dot.vector4_f32[2]
        vTemp = XM_PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
        // Result.vector4_f32[0] = (x+y)+z
        vDot = _mm_add_ss(vDot, vTemp);
        // Splat x
        return XM_PERMUTE_PS(vDot, _MM_SHUFFLE(0, 0, 0, 0));
#endif
}
    inline XMVECTOR XM_CALLCONV XMVector3LengthSq(FXMVECTOR V) noexcept
    {
        return XMVector3Dot(V, V);
    }
    inline XMVECTOR XM_CALLCONV XMVector3Length(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)

        XMVECTOR Result;

        Result = XMVector3LengthSq(V);
        Result = XMVectorSqrt(Result);

        return Result;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        // Dot3
        float32x4_t vTemp = vmulq_f32(V, V);
        float32x2_t v1 = vget_low_f32(vTemp);
        float32x2_t v2 = vget_high_f32(vTemp);
        v1 = vpadd_f32(v1, v1);
        v2 = vdup_lane_f32(v2, 0);
        v1 = vadd_f32(v1, v2);
        const float32x2_t zero = vdup_n_f32(0);
        uint32x2_t VEqualsZero = vceq_f32(v1, zero);
        // Sqrt
        float32x2_t S0 = vrsqrte_f32(v1);
        float32x2_t P0 = vmul_f32(v1, S0);
        float32x2_t R0 = vrsqrts_f32(P0, S0);
        float32x2_t S1 = vmul_f32(S0, R0);
        float32x2_t P1 = vmul_f32(v1, S1);
        float32x2_t R1 = vrsqrts_f32(P1, S1);
        float32x2_t Result = vmul_f32(S1, R1);
        Result = vmul_f32(v1, Result);
        Result = vbsl_f32(VEqualsZero, zero, Result);
        return vcombine_f32(Result, Result);
#elif defined(_XM_SSE4_INTRINSICS_)
        XMVECTOR vTemp = _mm_dp_ps(V, V, 0x7f);
        return _mm_sqrt_ps(vTemp);
#elif defined(_XM_SSE3_INTRINSICS_)
        XMVECTOR vLengthSq = _mm_mul_ps(V, V);
        vLengthSq = _mm_and_ps(vLengthSq, g_XMMask3);
        vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
        vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
        vLengthSq = _mm_sqrt_ps(vLengthSq);
        return vLengthSq;
#elif defined(_XM_SSE_INTRINSICS_)
        // Perform the dot product on x,y and z
        XMVECTOR vLengthSq = _mm_mul_ps(V, V);
        // vTemp has z and y
        XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(1, 2, 1, 2));
        // x+z, y
        vLengthSq = _mm_add_ss(vLengthSq, vTemp);
        // y,y,y,y
        vTemp = XM_PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
        // x+z+y,??,??,??
        vLengthSq = _mm_add_ss(vLengthSq, vTemp);
        // Splat the length squared
        vLengthSq = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
        // Get the length
        vLengthSq = _mm_sqrt_ps(vLengthSq);
        return vLengthSq;
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVector3Normalize(FXMVECTOR V) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        float fLength;
        XMVECTOR vResult;

        vResult = XMVector3Length(V);
        fLength = vResult.vector4_f32[0];

        // Prevent divide by zero
        if (fLength > 0)
        {
            fLength = 1.0f / fLength;
        }

        vResult.vector4_f32[0] = V.vector4_f32[0] * fLength;
        vResult.vector4_f32[1] = V.vector4_f32[1] * fLength;
        vResult.vector4_f32[2] = V.vector4_f32[2] * fLength;
        vResult.vector4_f32[3] = V.vector4_f32[3] * fLength;
        return vResult;

#elif defined(_XM_ARM_NEON_INTRINSICS_)
        // Dot3
        float32x4_t vTemp = vmulq_f32(V, V);
        float32x2_t v1 = vget_low_f32(vTemp);
        float32x2_t v2 = vget_high_f32(vTemp);
        v1 = vpadd_f32(v1, v1);
        v2 = vdup_lane_f32(v2, 0);
        v1 = vadd_f32(v1, v2);
        uint32x2_t VEqualsZero = vceq_f32(v1, vdup_n_f32(0));
        uint32x2_t VEqualsInf = vceq_f32(v1, vget_low_f32(g_XMInfinity));
        // Reciprocal sqrt (2 iterations of Newton-Raphson)
        float32x2_t S0 = vrsqrte_f32(v1);
        float32x2_t P0 = vmul_f32(v1, S0);
        float32x2_t R0 = vrsqrts_f32(P0, S0);
        float32x2_t S1 = vmul_f32(S0, R0);
        float32x2_t P1 = vmul_f32(v1, S1);
        float32x2_t R1 = vrsqrts_f32(P1, S1);
        v2 = vmul_f32(S1, R1);
        // Normalize
        XMVECTOR vResult = vmulq_f32(V, vcombine_f32(v2, v2));
        vResult = vbslq_f32(vcombine_f32(VEqualsZero, VEqualsZero), vdupq_n_f32(0), vResult);
        return vbslq_f32(vcombine_f32(VEqualsInf, VEqualsInf), g_XMQNaN, vResult);
#elif defined(_XM_SSE4_INTRINSICS_)
        XMVECTOR vLengthSq = _mm_dp_ps(V, V, 0x7f);
        // Prepare for the division
        XMVECTOR vResult = _mm_sqrt_ps(vLengthSq);
        // Create zero with a single instruction
        XMVECTOR vZeroMask = _mm_setzero_ps();
        // Test for a divide by zero (Must be FP to detect -0.0)
        vZeroMask = _mm_cmpneq_ps(vZeroMask, vResult);
        // Failsafe on zero (Or epsilon) length planes
        // If the length is infinity, set the elements to zero
        vLengthSq = _mm_cmpneq_ps(vLengthSq, g_XMInfinity);
        // Divide to perform the normalization
        vResult = _mm_div_ps(V, vResult);
        // Any that are infinity, set to zero
        vResult = _mm_and_ps(vResult, vZeroMask);
        // Select qnan or result based on infinite length
        XMVECTOR vTemp1 = _mm_andnot_ps(vLengthSq, g_XMQNaN);
        XMVECTOR vTemp2 = _mm_and_ps(vResult, vLengthSq);
        vResult = _mm_or_ps(vTemp1, vTemp2);
        return vResult;
#elif defined(_XM_SSE3_INTRINSICS_)
        // Perform the dot product on x,y and z only
        XMVECTOR vLengthSq = _mm_mul_ps(V, V);
        vLengthSq = _mm_and_ps(vLengthSq, g_XMMask3);
        vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
        vLengthSq = _mm_hadd_ps(vLengthSq, vLengthSq);
        // Prepare for the division
        XMVECTOR vResult = _mm_sqrt_ps(vLengthSq);
        // Create zero with a single instruction
        XMVECTOR vZeroMask = _mm_setzero_ps();
        // Test for a divide by zero (Must be FP to detect -0.0)
        vZeroMask = _mm_cmpneq_ps(vZeroMask, vResult);
        // Failsafe on zero (Or epsilon) length planes
        // If the length is infinity, set the elements to zero
        vLengthSq = _mm_cmpneq_ps(vLengthSq, g_XMInfinity);
        // Divide to perform the normalization
        vResult = _mm_div_ps(V, vResult);
        // Any that are infinity, set to zero
        vResult = _mm_and_ps(vResult, vZeroMask);
        // Select qnan or result based on infinite length
        XMVECTOR vTemp1 = _mm_andnot_ps(vLengthSq, g_XMQNaN);
        XMVECTOR vTemp2 = _mm_and_ps(vResult, vLengthSq);
        vResult = _mm_or_ps(vTemp1, vTemp2);
        return vResult;
#elif defined(_XM_SSE_INTRINSICS_)
        // Perform the dot product on x,y and z only
        XMVECTOR vLengthSq = _mm_mul_ps(V, V);
        XMVECTOR vTemp = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(2, 1, 2, 1));
        vLengthSq = _mm_add_ss(vLengthSq, vTemp);
        vTemp = XM_PERMUTE_PS(vTemp, _MM_SHUFFLE(1, 1, 1, 1));
        vLengthSq = _mm_add_ss(vLengthSq, vTemp);
        vLengthSq = XM_PERMUTE_PS(vLengthSq, _MM_SHUFFLE(0, 0, 0, 0));
        // Prepare for the division
        XMVECTOR vResult = _mm_sqrt_ps(vLengthSq);
        // Create zero with a single instruction
        XMVECTOR vZeroMask = _mm_setzero_ps();
        // Test for a divide by zero (Must be FP to detect -0.0)
        vZeroMask = _mm_cmpneq_ps(vZeroMask, vResult);
        // Failsafe on zero (Or epsilon) length planes
        // If the length is infinity, set the elements to zero
        vLengthSq = _mm_cmpneq_ps(vLengthSq, g_XMInfinity);
        // Divide to perform the normalization
        vResult = _mm_div_ps(V, vResult);
        // Any that are infinity, set to zero
        vResult = _mm_and_ps(vResult, vZeroMask);
        // Select qnan or result based on infinite length
        XMVECTOR vTemp1 = _mm_andnot_ps(vLengthSq, g_XMQNaN);
        XMVECTOR vTemp2 = _mm_and_ps(vResult, vLengthSq);
        vResult = _mm_or_ps(vTemp1, vTemp2);
        return vResult;
#endif
    }
    inline XMVECTOR XM_CALLCONV XMVector3Reflect(FXMVECTOR Incident, FXMVECTOR Normal) noexcept
    {
        // Result = Incident - (2 * dot(Incident, Normal)) * Normal

        XMVECTOR Result = XMVector3Dot(Incident, Normal);
        Result = XMVectorAdd(Result, Result);
        Result = XMVectorNegativeMultiplySubtract(Result, Normal, Incident);

        return Result;
}

    //--------------------Matrix------------------------------------------------
    inline XMMATRIX XM_CALLCONV XMMatrixMultiply(FXMMATRIX M1, CXMMATRIX M2) noexcept
    {
#if defined(_XM_NO_INTRINSICS_)
        XMMATRIX mResult;
        // Cache the invariants in registers
        float x = M1.m[0][0];
        float y = M1.m[0][1];
        float z = M1.m[0][2];
        float w = M1.m[0][3];
        // Perform the operation on the first row
        mResult.m[0][0] = (M2.m[0][0] * x) + (M2.m[1][0] * y) + (M2.m[2][0] * z) + (M2.m[3][0] * w);
        mResult.m[0][1] = (M2.m[0][1] * x) + (M2.m[1][1] * y) + (M2.m[2][1] * z) + (M2.m[3][1] * w);
        mResult.m[0][2] = (M2.m[0][2] * x) + (M2.m[1][2] * y) + (M2.m[2][2] * z) + (M2.m[3][2] * w);
        mResult.m[0][3] = (M2.m[0][3] * x) + (M2.m[1][3] * y) + (M2.m[2][3] * z) + (M2.m[3][3] * w);
        // Repeat for all the other rows
        x = M1.m[1][0];
        y = M1.m[1][1];
        z = M1.m[1][2];
        w = M1.m[1][3];
        mResult.m[1][0] = (M2.m[0][0] * x) + (M2.m[1][0] * y) + (M2.m[2][0] * z) + (M2.m[3][0] * w);
        mResult.m[1][1] = (M2.m[0][1] * x) + (M2.m[1][1] * y) + (M2.m[2][1] * z) + (M2.m[3][1] * w);
        mResult.m[1][2] = (M2.m[0][2] * x) + (M2.m[1][2] * y) + (M2.m[2][2] * z) + (M2.m[3][2] * w);
        mResult.m[1][3] = (M2.m[0][3] * x) + (M2.m[1][3] * y) + (M2.m[2][3] * z) + (M2.m[3][3] * w);
        x = M1.m[2][0];
        y = M1.m[2][1];
        z = M1.m[2][2];
        w = M1.m[2][3];
        mResult.m[2][0] = (M2.m[0][0] * x) + (M2.m[1][0] * y) + (M2.m[2][0] * z) + (M2.m[3][0] * w);
        mResult.m[2][1] = (M2.m[0][1] * x) + (M2.m[1][1] * y) + (M2.m[2][1] * z) + (M2.m[3][1] * w);
        mResult.m[2][2] = (M2.m[0][2] * x) + (M2.m[1][2] * y) + (M2.m[2][2] * z) + (M2.m[3][2] * w);
        mResult.m[2][3] = (M2.m[0][3] * x) + (M2.m[1][3] * y) + (M2.m[2][3] * z) + (M2.m[3][3] * w);
        x = M1.m[3][0];
        y = M1.m[3][1];
        z = M1.m[3][2];
        w = M1.m[3][3];
        mResult.m[3][0] = (M2.m[0][0] * x) + (M2.m[1][0] * y) + (M2.m[2][0] * z) + (M2.m[3][0] * w);
        mResult.m[3][1] = (M2.m[0][1] * x) + (M2.m[1][1] * y) + (M2.m[2][1] * z) + (M2.m[3][1] * w);
        mResult.m[3][2] = (M2.m[0][2] * x) + (M2.m[1][2] * y) + (M2.m[2][2] * z) + (M2.m[3][2] * w);
        mResult.m[3][3] = (M2.m[0][3] * x) + (M2.m[1][3] * y) + (M2.m[2][3] * z) + (M2.m[3][3] * w);
        return mResult;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
        XMMATRIX mResult;
        float32x2_t VL = vget_low_f32(M1.r[0]);
        float32x2_t VH = vget_high_f32(M1.r[0]);
        // Perform the operation on the first row
        XMVECTOR vX = vmulq_lane_f32(M2.r[0], VL, 0);
        XMVECTOR vY = vmulq_lane_f32(M2.r[1], VL, 1);
        XMVECTOR vZ = vmlaq_lane_f32(vX, M2.r[2], VH, 0);
        XMVECTOR vW = vmlaq_lane_f32(vY, M2.r[3], VH, 1);
        mResult.r[0] = vaddq_f32(vZ, vW);
        // Repeat for the other 3 rows
        VL = vget_low_f32(M1.r[1]);
        VH = vget_high_f32(M1.r[1]);
        vX = vmulq_lane_f32(M2.r[0], VL, 0);
        vY = vmulq_lane_f32(M2.r[1], VL, 1);
        vZ = vmlaq_lane_f32(vX, M2.r[2], VH, 0);
        vW = vmlaq_lane_f32(vY, M2.r[3], VH, 1);
        mResult.r[1] = vaddq_f32(vZ, vW);
        VL = vget_low_f32(M1.r[2]);
        VH = vget_high_f32(M1.r[2]);
        vX = vmulq_lane_f32(M2.r[0], VL, 0);
        vY = vmulq_lane_f32(M2.r[1], VL, 1);
        vZ = vmlaq_lane_f32(vX, M2.r[2], VH, 0);
        vW = vmlaq_lane_f32(vY, M2.r[3], VH, 1);
        mResult.r[2] = vaddq_f32(vZ, vW);
        VL = vget_low_f32(M1.r[3]);
        VH = vget_high_f32(M1.r[3]);
        vX = vmulq_lane_f32(M2.r[0], VL, 0);
        vY = vmulq_lane_f32(M2.r[1], VL, 1);
        vZ = vmlaq_lane_f32(vX, M2.r[2], VH, 0);
        vW = vmlaq_lane_f32(vY, M2.r[3], VH, 1);
        mResult.r[3] = vaddq_f32(vZ, vW);
        return mResult;
#elif defined(_XM_AVX2_INTRINSICS_)
        __m256 t0 = _mm256_castps128_ps256(M1.r[0]);
        t0 = _mm256_insertf128_ps(t0, M1.r[1], 1);
        __m256 t1 = _mm256_castps128_ps256(M1.r[2]);
        t1 = _mm256_insertf128_ps(t1, M1.r[3], 1);

        __m256 u0 = _mm256_castps128_ps256(M2.r[0]);
        u0 = _mm256_insertf128_ps(u0, M2.r[1], 1);
        __m256 u1 = _mm256_castps128_ps256(M2.r[2]);
        u1 = _mm256_insertf128_ps(u1, M2.r[3], 1);

        __m256 a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(0, 0, 0, 0));
        __m256 a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(0, 0, 0, 0));
        __m256 b0 = _mm256_permute2f128_ps(u0, u0, 0x00);
        __m256 c0 = _mm256_mul_ps(a0, b0);
        __m256 c1 = _mm256_mul_ps(a1, b0);

        a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(1, 1, 1, 1));
        a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(1, 1, 1, 1));
        b0 = _mm256_permute2f128_ps(u0, u0, 0x11);
        __m256 c2 = _mm256_fmadd_ps(a0, b0, c0);
        __m256 c3 = _mm256_fmadd_ps(a1, b0, c1);

        a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(2, 2, 2, 2));
        a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(2, 2, 2, 2));
        __m256 b1 = _mm256_permute2f128_ps(u1, u1, 0x00);
        __m256 c4 = _mm256_mul_ps(a0, b1);
        __m256 c5 = _mm256_mul_ps(a1, b1);

        a0 = _mm256_shuffle_ps(t0, t0, _MM_SHUFFLE(3, 3, 3, 3));
        a1 = _mm256_shuffle_ps(t1, t1, _MM_SHUFFLE(3, 3, 3, 3));
        b1 = _mm256_permute2f128_ps(u1, u1, 0x11);
        __m256 c6 = _mm256_fmadd_ps(a0, b1, c4);
        __m256 c7 = _mm256_fmadd_ps(a1, b1, c5);

        t0 = _mm256_add_ps(c2, c6);
        t1 = _mm256_add_ps(c3, c7);

        XMMATRIX mResult;
        mResult.r[0] = _mm256_castps256_ps128(t0);
        mResult.r[1] = _mm256_extractf128_ps(t0, 1);
        mResult.r[2] = _mm256_castps256_ps128(t1);
        mResult.r[3] = _mm256_extractf128_ps(t1, 1);
        return mResult;
#elif defined(_XM_SSE_INTRINSICS_)
        XMMATRIX mResult;
        // Splat the component X,Y,Z then W
#if defined(_XM_AVX_INTRINSICS_)
        XMVECTOR vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 0);
        XMVECTOR vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 1);
        XMVECTOR vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 2);
        XMVECTOR vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[0]) + 3);
#else
    // Use vW to hold the original row
        XMVECTOR vW = M1.r[0];
        XMVECTOR vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
        XMVECTOR vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
        XMVECTOR vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
        vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
        // Perform the operation on the first row
        vX = _mm_mul_ps(vX, M2.r[0]);
        vY = _mm_mul_ps(vY, M2.r[1]);
        vZ = _mm_mul_ps(vZ, M2.r[2]);
        vW = _mm_mul_ps(vW, M2.r[3]);
        // Perform a binary add to reduce cumulative errors
        vX = _mm_add_ps(vX, vZ);
        vY = _mm_add_ps(vY, vW);
        vX = _mm_add_ps(vX, vY);
        mResult.r[0] = vX;
        // Repeat for the other 3 rows
#if defined(_XM_AVX_INTRINSICS_)
        vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 0);
        vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 1);
        vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 2);
        vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[1]) + 3);
#else
        vW = M1.r[1];
        vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
        vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
        vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
        vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
        vX = _mm_mul_ps(vX, M2.r[0]);
        vY = _mm_mul_ps(vY, M2.r[1]);
        vZ = _mm_mul_ps(vZ, M2.r[2]);
        vW = _mm_mul_ps(vW, M2.r[3]);
        vX = _mm_add_ps(vX, vZ);
        vY = _mm_add_ps(vY, vW);
        vX = _mm_add_ps(vX, vY);
        mResult.r[1] = vX;
#if defined(_XM_AVX_INTRINSICS_)
        vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 0);
        vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 1);
        vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 2);
        vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[2]) + 3);
#else
        vW = M1.r[2];
        vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
        vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
        vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
        vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
        vX = _mm_mul_ps(vX, M2.r[0]);
        vY = _mm_mul_ps(vY, M2.r[1]);
        vZ = _mm_mul_ps(vZ, M2.r[2]);
        vW = _mm_mul_ps(vW, M2.r[3]);
        vX = _mm_add_ps(vX, vZ);
        vY = _mm_add_ps(vY, vW);
        vX = _mm_add_ps(vX, vY);
        mResult.r[2] = vX;
#if defined(_XM_AVX_INTRINSICS_)
        vX = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 0);
        vY = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 1);
        vZ = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 2);
        vW = _mm_broadcast_ss(reinterpret_cast<const float*>(&M1.r[3]) + 3);
#else
        vW = M1.r[3];
        vX = XM_PERMUTE_PS(vW, _MM_SHUFFLE(0, 0, 0, 0));
        vY = XM_PERMUTE_PS(vW, _MM_SHUFFLE(1, 1, 1, 1));
        vZ = XM_PERMUTE_PS(vW, _MM_SHUFFLE(2, 2, 2, 2));
        vW = XM_PERMUTE_PS(vW, _MM_SHUFFLE(3, 3, 3, 3));
#endif
        vX = _mm_mul_ps(vX, M2.r[0]);
        vY = _mm_mul_ps(vY, M2.r[1]);
        vZ = _mm_mul_ps(vZ, M2.r[2]);
        vW = _mm_mul_ps(vW, M2.r[3]);
        vX = _mm_add_ps(vX, vZ);
        vY = _mm_add_ps(vY, vW);
        vX = _mm_add_ps(vX, vY);
        mResult.r[3] = vX;
        return mResult;
#endif
    }
    inline XMMATRIX XM_CALLCONV XMMATRIX::operator*(FXMMATRIX M) const noexcept
    {
        return XMMatrixMultiply(*this, M);
    }
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#pragma warning(pop)

    //------------------------------------------------------------------------------
    // Vector operators

#ifndef _XM_NO_XMVECTOR_OVERLOADS_
    inline XMVECTOR    XM_CALLCONV     operator+ (FXMVECTOR V) noexcept
    {
        return V;
    }
    inline XMVECTOR    XM_CALLCONV     operator- (FXMVECTOR V) noexcept
    {
        return XMVectorNegate(V);
    }

    inline XMVECTOR& XM_CALLCONV     operator+= (XMVECTOR& V1, FXMVECTOR V2) noexcept
    {
        V1 = XMVectorAdd(V1, V2);
        return V1;
    }
    inline XMVECTOR& XM_CALLCONV     operator*= (XMVECTOR& V1, FXMVECTOR V2) noexcept
    {
        V1 = XMVectorMultiply(V1, V2);
        return V1;
    }

    inline XMVECTOR    XM_CALLCONV     operator+ (FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
        return XMVectorAdd(V1, V2);
    }
    inline XMVECTOR    XM_CALLCONV     operator- (FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
        return XMVectorSubtract(V1, V2);
    }
    inline XMVECTOR    XM_CALLCONV     operator* (FXMVECTOR V1, FXMVECTOR V2) noexcept
    {
        return XMVectorMultiply(V1, V2);
    }
#endif /* !_XM_NO_XMVECTOR_OVERLOADS_ */

    namespace PackedVector
    {
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-anonymous-struct"
#pragma clang diagnostic ignored "-Wnested-anon-types"
#endif
        struct XMCOLOR
        {
            union
            {
                struct
                {
                    uint8_t b;  // Blue:    0/255 to 255/255
                    uint8_t g;  // Green:   0/255 to 255/255
                    uint8_t r;  // Red:     0/255 to 255/255
                    uint8_t a;  // Alpha:   0/255 to 255/255
                };
                uint32_t c;
            };
        public:
            XMCOLOR() = default;

            XMCOLOR(const XMCOLOR&) = default;
            XMCOLOR& operator=(const XMCOLOR&) = default;

            XMCOLOR(XMCOLOR&&) = default;
            XMCOLOR& operator=(XMCOLOR&&) = default;
            constexpr XMCOLOR(uint32_t Color) noexcept : c(Color) {}
        public:
            operator uint32_t () const noexcept { return c; }
            XMCOLOR& operator= (const uint32_t Color) noexcept { c = Color; return *this; }
        };

 
#ifdef __clang__
#pragma clang diagnostic pop
#endif
        inline void XM_CALLCONV XMStoreColor2
        (
            DirectX::PackedVector::XMCOLOR* pDestination,
            DirectX::FXMVECTOR V
        ) noexcept
        {
            using namespace DirectX;
            assert(pDestination);
#if defined(_XM_NO_INTRINSICS_)

            XMVECTOR N = XMVectorSaturate(V);
            N = XMVectorMultiply(N, g_UByteMax);
            N = XMVectorRound(N);

            XMFLOAT4A tmp;
            XMStoreFloat4A(&tmp, N);

            pDestination->c = (static_cast<uint32_t>(tmp.w) << 24) |
                (static_cast<uint32_t>(tmp.x) << 16) |
                (static_cast<uint32_t>(tmp.y) << 8) |
                static_cast<uint32_t>(tmp.z);

#elif defined(_XM_ARM_NEON_INTRINSICS_)
            float32x4_t R = vmaxq_f32(V, vdupq_n_f32(0));
            R = vminq_f32(R, vdupq_n_f32(1.0f));
            R = vmulq_n_f32(R, 255.0f);
            R = XMVectorRound(R);
            uint32x4_t vInt32 = vcvtq_u32_f32(R);
            uint16x4_t vInt16 = vqmovn_u32(vInt32);
            uint8x8_t vInt8 = vqmovn_u16(vcombine_u16(vInt16, vInt16));
            uint32_t rgba = vget_lane_u32(vreinterpret_u32_u8(vInt8), 0);
            pDestination->c = (rgba & 0xFF00FF00) | ((rgba >> 16) & 0xFF) | ((rgba << 16) & 0xFF0000);
#elif defined(_XM_SSE_INTRINSICS_)
            // Set <0 to 0
            XMVECTOR vResult = _mm_max_ps(V, g_XMZero);
            // Set>1 to 1
            vResult = _mm_min_ps(vResult, g_XMOne);
             //Convert to 0-255
            vResult = _mm_mul_ps(vResult, g_UByteMax);
             //Convert to int
            __m128i vInt = _mm_cvtps_epi32(vResult);
            // Mash to shorts
            vInt = _mm_packs_epi32(vInt, vInt);
            // Mash to bytes
            vInt = _mm_packus_epi16(vInt, vInt);
            // Store the color
            _mm_store_ss(reinterpret_cast<float*>(&pDestination->c), _mm_castsi128_ps(vInt));
#endif
        }
    } // namespace PackedVector

} // namespace DirectX

#pragma endregion
#pragma region XMVS
union alignas(16)VMFLOAT32A
{
    VMFLOAT32A() = default;
    VMFLOAT32A(__m128 in) :v(in) {};
public:
    VMFLOAT32A& operator+=(const VMFLOAT32A & in)
    {
        using namespace DirectX;
        v = (v + in.v);
        return *this;
    }
    operator DirectX::XMVECTOR()const
    {
        return v;
    }
    VMFLOAT32A& operator=(const __m128 in)
    {
        v = in;
        return *this;
    }
public:
    float f[4];
    DirectX::XMVECTOR v;
};
union alignas(16)XMVSOut
{
public:
    XMVSOut() :gl{} {};
public:
    void Increase(const XMVSOut & other, size_t size)
    {
        for (size_t i = 0; i < size; i++)
        {
            dx.attributes[i] += other.dx.attributes[i];
        }
        dx.SV_Position += other.dx.SV_Position;
    }
    void Scale(DirectX::XMVECTOR by, size_t size)
    {
        using namespace DirectX;
        dx.SV_Position.v *= by;
        for (size_t i = 0; i < size; i++)
        {
            dx.attributes[i].v *= by;
        }
    }

    XMVSOut& operator+=(const XMVSOut & rhs)
    {
        using namespace DirectX;
        dx.SV_Position.v += rhs.dx.SV_Position.v;
        for (size_t i = 0; i < maxAttributes; i++)
        {
            dx.attributes[i].v += rhs.dx.attributes[i].v;
        }
        return *this;
    }
    static XMVSOut Subtract(const XMVSOut & lhs, const XMVSOut & rhs, size_t size)
    {
        using namespace DirectX;
        XMVSOut out;
        out.dx.SV_Position.v = lhs.dx.SV_Position.v - rhs.dx.SV_Position.v;
        for (size_t i = 0; i < size; i++)
        {
            out.dx.attributes[i].v = lhs.dx.attributes[i].v - rhs.dx.attributes[i].v;
        }
        return out;
    }
    static XMVSOut Multiply(const XMVSOut & lhs, DirectX::FXMVECTOR rhs, size_t size)
    {
        using namespace DirectX;
        XMVSOut out;
        out.dx.SV_Position.v = lhs.dx.SV_Position.v * rhs;
        for (size_t i = 0; i < size; i++)
        {
            out.dx.attributes[i].v = lhs.dx.attributes[i].v * rhs;
        }
        return out;
    }
public:
    struct DXVSOut
    {
        VMFLOAT32A attributes[maxAttributes];
        VMFLOAT32A SV_Position;
    }dx;
    OutVertex gl;
};
union alignas(16)XMOutPixel
{
    XMOutPixel() : SV_Target{} {};
    DirectX::XMVECTORF32 SV_Target;
    OutFragment gl;
};
union alignas(16)XMInPixel
{
    XMInPixel() :dx{} {};
    XMVSOut dx;
    InFragment gl;
};
typedef const VMFLOAT32A FVMFLOAT32A;


inline XMVSOut VSOutInterpolate(const XMVSOut& v0, const XMVSOut& v1, float alpha, size_t voSize)
{
    using namespace DirectX;
    XMVSOut out;
    for (size_t i = 0; i < voSize; i++)
    {
        out.dx.attributes[i].v = XMVectorLerp(v0.dx.attributes[i], v1.dx.attributes[i], alpha);
    }
    out.dx.SV_Position.v = XMVectorLerp(v0.dx.SV_Position, v1.dx.SV_Position, alpha);
    return out;
}

inline VMFLOAT32A XM_CALLCONV VMVector3Cross
(
    FVMFLOAT32A V1,
    FVMFLOAT32A V2
) noexcept
{
    // [ V1.y*V2.z - V1.z*V2.y, V1.z*V2.x - V1.x*V2.z, V1.x*V2.y - V1.y*V2.x ]

#if defined(_XM_NO_INTRINSICS_)
    XMVECTORF32 vResult = { { {
            (V1.vector4_f32[1] * V2.vector4_f32[2]) - (V1.vector4_f32[2] * V2.vector4_f32[1]),
            (V1.vector4_f32[2] * V2.vector4_f32[0]) - (V1.vector4_f32[0] * V2.vector4_f32[2]),
            (V1.vector4_f32[0] * V2.vector4_f32[1]) - (V1.vector4_f32[1] * V2.vector4_f32[0]),
            0.0f
        } } };
    return vResult.v;
#elif defined(_XM_ARM_NEON_INTRINSICS_)
    float32x2_t v1xy = vget_low_f32(V1);
    float32x2_t v2xy = vget_low_f32(V2);

    float32x2_t v1yx = vrev64_f32(v1xy);
    float32x2_t v2yx = vrev64_f32(v2xy);

    float32x2_t v1zz = vdup_lane_f32(vget_high_f32(V1), 0);
    float32x2_t v2zz = vdup_lane_f32(vget_high_f32(V2), 0);

    XMVECTOR vResult = vmulq_f32(vcombine_f32(v1yx, v1xy), vcombine_f32(v2zz, v2yx));
    vResult = vmlsq_f32(vResult, vcombine_f32(v1zz, v1yx), vcombine_f32(v2yx, v2xy));
    vResult = veorq_u32(vResult, g_XMFlipY);
    return vandq_u32(vResult, g_XMMask3);
#elif defined(_XM_SSE_INTRINSICS_)
    // y1,z1,x1,w1
    DirectX::XMVECTOR vTemp1 = XM_PERMUTE_PS(V1, _MM_SHUFFLE(3, 0, 2, 1));
    // z2,x2,y2,w2
    DirectX::XMVECTOR vTemp2 = XM_PERMUTE_PS(V2, _MM_SHUFFLE(3, 1, 0, 2));
    // Perform the left operation
    DirectX::XMVECTOR vResult = _mm_mul_ps(vTemp1, vTemp2);
    // z1,x1,y1,w1
    vTemp1 = XM_PERMUTE_PS(vTemp1, _MM_SHUFFLE(3, 0, 2, 1));
    // y2,z2,x2,w2
    vTemp2 = XM_PERMUTE_PS(vTemp2, _MM_SHUFFLE(3, 1, 0, 2));
    // Perform the right operation
    vResult = XM_FNMADD_PS(vTemp1, vTemp2, vResult);
    // Set w to zero
    return _mm_and_ps(vResult, DirectX::g_XMMask3);
#endif
}
#pragma endregion


namespace vc
{
    class DragRace
    {
    public:
        DragRace(size_t thread_count = std::thread::hardware_concurrency())
        {
            threads.reserve(thread_count);
            tasks.resize(thread_count);
            results.resize(thread_count);


            for (size_t i = 0; i < thread_count; i++)
            {
                threads.emplace_back(
                    [&]() 
                    {
                        std::unique_lock<std::mutex> doorlock(task_mutex, std::defer_lock);

                        while (true)
                        {
                            doorlock.lock();
                            task_cv.wait(doorlock,
                                [&]() -> bool { return tasks_set != 0 || stop_threads; }
                            );
                            if (stop_threads && !tasks_set) return;
                            auto index = --tasks_set; 
                            doorlock.unlock();

                            
                            auto temp_task = std::move(tasks[index]);
                            temp_task();
                        }
                    }
                );
            }
        }
        ~DragRace() 
        {
            stop_threads = true;
            task_cv.notify_all();

            for (std::thread& thread : threads)
            {
                thread.join();
            }
        }
        DragRace(const DragRace&) = delete;
        DragRace& operator=(const DragRace&) = delete;
    public:
        template <typename F>
        void LoadTask(F&& function)
        {
            std::packaged_task<std::invoke_result_t<F>()> task_pkg(std::forward<F>(function));
            results[tasks_set] = task_pkg.get_future();
            tasks[tasks_set++] = std::packaged_task<void()>(std::move(task_pkg));
            if (tasks_set == (uint8_t)ThreadCount()) Launch();
        }
        size_t ThreadCount() const noexcept
        {
            return threads.size();
        }
    private:
        void Launch()
        {
            task_cv.notify_all();
            for (auto& th : results)
            {
                th.wait();
            }
        }
    private:
        std::vector<std::thread> threads; //racers(horses)
        std::vector<std::packaged_task<void()>> tasks; //stakes
        std::vector<std::future<void>> results;

        std::mutex task_mutex;
        std::condition_variable task_cv;

        bool stop_threads = false;
        uint8_t tasks_set = 0;
    };
}
class Buffer
{
public:
    Buffer(size_t size)
    {
        Body.resize(size);
    }
    Buffer(Buffer&& in)noexcept
        :Body(std::move(in.Body))
    {
    }
public:
    BufferID GetID()const noexcept
    {
        return reinterpret_cast<BufferID>(&Body[0]);
    }
    uint8_t* GetRawData()noexcept
    {
        return &Body[0];
    }
    void InsertDataChunk(const void* data, size_t offset, size_t volume)
    {
        std::copy(static_cast<const uint8_t*>(data), static_cast<const uint8_t*>(data) + volume, Body.begin() + offset);
    }
    template<typename Type>
    constexpr const std::vector<Type>& Cast()const
    {
        return reinterpret_cast<const std::vector<Type>&>(Body);
    }
    template<typename T>
    constexpr const T& Get(size_t n, size_t stride, size_t offset)const
    {
        return reinterpret_cast<const T&>(Body[n * stride + offset]);
    }
private:
    std::vector<uint8_t> Body;
};
class InputLayout
{
public:
    struct InputLayoutDesc
    {
        InputLayoutDesc() = default;
        InputLayoutDesc(const Buffer* buff, size_t Offset, size_t Stride, AttributeType Format)
            :CurBuffer(buff), Offset(Offset), Stride(Stride), Format(Format)
        {}
    public:
        const Buffer* CurBuffer = nullptr;
        size_t Offset = 0;
        size_t Stride = 0;
        AttributeType Format = AttributeType::EMPTY;
    };
public:
    InputLayout() {};
public:
    void EnableHead(size_t n)
    {
        ActiveHeads[n] = true;
    }
    void DisableHead(size_t n)
    {
        ActiveHeads[n] = false;
    }
    void SetIndexing(const Buffer* IndexBuffer, IndexType Mode = IndexType::UINT32)
    {
        IAActiveIndexBuffer = IndexBuffer;
        indexing = Mode;
    }
    void SetDescriptor(size_t head, AttributeType type, size_t stride, size_t offset, const Buffer* buffer)
    {
        descriptors[head] = { buffer, offset, stride, type };
    }
    bool IsIndexed()const
    {
        return IAActiveIndexBuffer != nullptr;
    }
    InVertex MakeVertex(size_t SV_VertexID)const
    {
        InVertex VSInVert;
        VSInVert.gl_VertexID = (uint32_t)SV_VertexID;
        for (auto i = 0u; i < maxAttributes; i++)
        {
            if (ActiveHeads[i])
            {
                auto& desc = descriptors[i];
                switch (desc.Format)
                {
                case AttributeType::FLOAT:
                {
                    VSInVert.attributes[i].v1 = desc.CurBuffer->Get<float>(SV_VertexID, desc.Stride, desc.Offset); break;
                }
                case AttributeType::VEC2:
                {
                    VSInVert.attributes[i].v2 = desc.CurBuffer->Get<glm::vec2>(SV_VertexID, desc.Stride, desc.Offset); break;
                }
                case AttributeType::VEC3:
                {
                    VSInVert.attributes[i].v3 = desc.CurBuffer->Get<glm::vec3>(SV_VertexID, desc.Stride, desc.Offset); break;
                }
                case AttributeType::VEC4:
                {
                    VSInVert.attributes[i].v4 = desc.CurBuffer->Get<glm::vec4>(SV_VertexID, desc.Stride, desc.Offset); break;
                }
                default:
                {
                    assert(false && "Head was enabled, but no matching type was found");
                    break;
                }
                }
            }
        }
        return VSInVert;
    }
    template<typename T>
    const auto& GetIndexBuffer()const
    {
        return IAActiveIndexBuffer->Cast<T>();
    }
    IndexType GetIndexType()const noexcept
    {
        return indexing;
    }
private:
    std::array<InputLayoutDesc, maxAttributes> descriptors;
    const Buffer* IAActiveIndexBuffer = nullptr;
    IndexType indexing = IndexType::UINT32;
    std::bitset<maxAttributes> ActiveHeads;
};
class ShaderProgramm
{
    friend class GPU;
    union Constant
    {
        Constant() {};
        Constant(float v1) :v1{ v1 } {};
        Constant(glm::vec2 v2) :v2{ v2 } {};
        Constant(glm::vec3 v3) :v3{ v3 } {};
        Constant(glm::vec4 v4) :v4{ v4 } {};
        Constant(glm::mat4x4 m4) :m4{ m4 } {};
    public:
        float     v1; ///< single float
        glm::vec2 v2; ///< two floats
        glm::vec3 v3; ///< three floats
        glm::vec4 v4; ///< four floats
        glm::mat4 m4 = glm::mat4(1.f); ///< 4x4 float matrix
    };
    using VertexShaderEX = std::function<void(OutVertex&, const InVertex&, const Uniforms&)>;
    using PixelShaderEX = std::function<void(OutFragment&, const InFragment&, const Uniforms&)>;
    using ConstantBuffer = std::array<Constant, maxAttributes>;
public:
    ShaderProgramm() = default;
    ShaderProgramm(PixelShaderEX PS_in, VertexShaderEX VS_in) :PS(PS_in), VS(VS_in) {}
public:
    void SetPixelShader(PixelShaderEX PS_in)noexcept
    {
        PS = PS_in;
    }
    void SetVertexShader(VertexShaderEX VS_in)noexcept
    {
        VS = VS_in;
    }
    void SetVStoPSAttibutes(size_t index, AttributeType attr)noexcept
    {
        VStoPS[index] = attr;
        if (MonotonicSize < index + 1)
        {
            MonotonicSize = index + 1;
        }
    }
    void SetVStoPSAttibutes(std::initializer_list<AttributeType> attrs, size_t start_index = 0)noexcept
    {
        for (auto attr : attrs)
        {
            VStoPS[start_index] = attr;
            start_index++;
        }
        MonotonicSize = start_index;
    }
    size_t GetMonotonicSize()const noexcept
    {
        return MonotonicSize;
    }

    XMOutPixel InvokePS(const InFragment& in)
    {
        XMOutPixel _out;
        PS(_out.gl, in, (Uniforms&)ConstantBuffers[0]);
        return _out;
    }
    void InvokeVS(OutVertex& out, const InVertex& in)
    {
        VS(out, in, (Uniforms&)ConstantBuffers[0]);
    }

    template<typename MathT> constexpr
        void SetConstantBuffer(size_t index, MathT&& c)
    {
        ConstantBuffers[index] = c;
    }
    template<typename ...MathT> constexpr
        void SetConstantBuffer(size_t offset, MathT&&...c)
    {
        (SetConstantBuffer(offset++, std::forward<MathT&&>(c)), ...);
    }
private:
    PixelShaderEX PS;
    VertexShaderEX VS;
    ConstantBuffer ConstantBuffers;
    std::array<AttributeType, maxAttributes> VStoPS{ AttributeType::EMPTY };
    size_t MonotonicSize = 0u;
};

class GPU final
{
public:
    GPU();
public:
    BufferID  createBuffer           (uint64_t size);
    BufferID  createBuffer           (Buffer** _out, uint64_t size);
    void      deleteBuffer           (BufferID buffer);
    void      setBufferData          (BufferID buffer, uint64_t offset, uint64_t size, const void* data);
    void      getBufferData          (BufferID buffer, uint64_t offset, uint64_t size, void* data);
    bool      isBuffer               (BufferID buffer);

public:    //vertex array object commands (vertex puller)
    ObjectID  createVertexPuller     ();
    ObjectID  createVertexPuller     (InputLayout** _out);
    void      deleteVertexPuller     (VertexPullerID vao);
    void      setVertexPullerHead    (VertexPullerID vao,uint32_t head,AttributeType type,uint64_t stride,uint64_t offset,BufferID buffer);
    void      setVertexPullerIndexing(VertexPullerID vao,IndexType type,BufferID buffer);
    void      enableVertexPullerHead (VertexPullerID vao,uint32_t head);
    void      disableVertexPullerHead(VertexPullerID vao,uint32_t head);
    void      bindVertexPuller       (VertexPullerID vao);
    void      bindVertexPuller       (InputLayout* in);
    void      unbindVertexPuller     ();
    bool      isVertexPuller         (VertexPullerID vao);

    //program object commands
    ProgramID createProgram          ();
    ProgramID createProgram          (ShaderProgramm** _out, VertexShader vs, FragmentShader fs);
    void      deleteProgram          (ProgramID prg);
    void      attachShaders          (ProgramID prg,VertexShader vs,FragmentShader fs);
    void      setVS2FSType           (ProgramID prg,uint32_t attrib,AttributeType type);
    void      useProgram             (ProgramID prg);
    void      useProgram             (ShaderProgramm* prg);
    bool      isProgram              (ProgramID prg);
    void      programUniform1f       (ProgramID prg,uint32_t uniformId,float     const&d);
    void      programUniform2f       (ProgramID prg,uint32_t uniformId,glm::vec2 const&d);
    void      programUniform3f       (ProgramID prg,uint32_t uniformId,glm::vec3 const&d);
    void      programUniform4f       (ProgramID prg,uint32_t uniformId,glm::vec4 const&d);
    void      programUniformMatrix4f (ProgramID prg,uint32_t uniformId,glm::mat4 const&d);

    //framebuffer functions
    void      createFramebuffer      (uint32_t width,uint32_t height);
    void      deleteFramebuffer      ();
    void      resizeFramebuffer      (uint32_t width,uint32_t height);
    uint8_t*  getFramebufferColor    ();
    float*    getFramebufferDepth    ();
    uint32_t  getFramebufferWidth    ();
    uint32_t  getFramebufferHeight   ();

    //execution commands
    void      clear                  (float r,float g,float b,float a);
    void      drawTriangles          (uint32_t  nofVertices);

    template<typename IndexType = uint32_t>
    void drawIndexedNoMT(uint32_t nofVertices)
    {
        assert(nofVertices % 3 == 0 && "Size is not a multiplier of 3");
        auto& Assembler = *activeIA;
        auto& Program = *activeSP;

        const auto& indexArray = Assembler.GetIndexBuffer<IndexType>();
        
       
        std::vector<XMVSOut> VSOut{ size_t(nofVertices) };
        InVertex VSInVert;
        uint32_t SV_VertexID = 0;

        for (size_t i = 0; i < nofVertices; i += 3)
        {
            auto& v0 = VSOut[i];
            auto& v1 = VSOut[i + 1];
            auto& v2 = VSOut[i + 2];

            SV_VertexID = uint32_t(indexArray[i]);
            Program.InvokeVS(v0.gl, Assembler.MakeVertex(SV_VertexID));
            SV_VertexID = uint32_t(indexArray[i + 1]);
            Program.InvokeVS(v1.gl, Assembler.MakeVertex(SV_VertexID));
            SV_VertexID = uint32_t(indexArray[i + 2]);
            Program.InvokeVS(v2.gl, Assembler.MakeVertex(SV_VertexID));

            ClipCullTriangles(v0, v1, v2);
        };
    }
    template<typename IndexType = uint32_t>
    void drawIndexed(uint32_t nofVertices)
    {
        assert(nofVertices % 3 == 0 && "Size is not a multiplier of 3");
        if (nofVertices < 3 * race.ThreadCount())
        {
            drawIndexedNoMT<IndexType>(nofVertices);
            return;
        }
        auto& Assembler = *activeIA;
        auto& Program = *activeSP;

        const auto& indexArray = Assembler.GetIndexBuffer<IndexType>();
        size_t begin = 0;
        size_t delta = nofVertices / race.ThreadCount() + 3 - (nofVertices / race.ThreadCount()) % 3;
        size_t end = delta;

        auto plF = [&](size_t begin, size_t end)
        {
            uint32_t SV_VertexID = 0;
            for (size_t i = begin; i < end; i += 3)
            {
                XMVSOut v0;
                XMVSOut v1;
                XMVSOut v2;

                SV_VertexID = uint32_t(indexArray[i]);
                Program.InvokeVS(v0.gl, Assembler.MakeVertex(SV_VertexID));
                SV_VertexID = uint32_t(indexArray[i + 1]);
                Program.InvokeVS(v1.gl, Assembler.MakeVertex(SV_VertexID));
                SV_VertexID = uint32_t(indexArray[i + 2]);
                Program.InvokeVS(v2.gl, Assembler.MakeVertex(SV_VertexID));

                ClipCullTriangles(v0, v1, v2);
            }
        };

        for (size_t i = 1; i < race.ThreadCount(); i++)
        {
            race.LoadTask(std::bind(plF, begin, end));
            begin = end;
            end += delta;
        }
        race.LoadTask(std::bind(plF, begin, nofVertices));
    }
private:
    void AssembleTriangles(std::vector<XMVSOut>& VSOut);
    void ClipCullTriangles(XMVSOut& Vo1, XMVSOut& Vo2, XMVSOut& Vo3);
    void PostProcessTriangle(XMVSOut& Vo1, XMVSOut& Vo2, XMVSOut& Vo3);
    void DrawTriangle(const XMVSOut& Vo1, const XMVSOut& Vo2, const XMVSOut& Vo3);
    void DrawFlatTopTriangle(const XMVSOut& Vo1, const XMVSOut& Vo2, const XMVSOut& Vo3);
    void DrawFlatBottomTriangle(const XMVSOut& Vo1, const XMVSOut& Vo2, const XMVSOut& Vo3);
    void DrawFlatTriangle(const XMVSOut& it0, const XMVSOut& it2, const XMVSOut& dv0, const XMVSOut& dv1, XMVSOut& itEdge1);
    std::pair<bool, float> DepthTest(uint32_t width_in, size_t PremulIndex, float z);
private:
    vc::DragRace race{ 16 };
    std::unordered_map<BufferID, Buffer> buffers;
    std::unordered_map<VertexPullerID, std::unique_ptr<InputLayout>> InputAssemblers;
    std::unordered_map<ProgramID, std::unique_ptr<ShaderProgramm>> ShaderStage;
    
    std::vector<DirectX::PackedVector::XMCOLOR> FrameBuffer;
    std::vector<float> DepthStencil;

    VMFLOAT32A Scale;
    VMFLOAT32A Offset;

    InputLayout* activeIA = nullptr;
    ShaderProgramm* activeSP = nullptr;

    uint32_t width = 0;
    uint32_t height = 0;
};




