//
// Copyright (C) YuqiaoZhang(HanetakaChou)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

#ifndef _BRX_FMINSEARCH_H_
#define _BRX_FMINSEARCH_H_ 1

static inline float internal_brx_fminsearch_3d_guarded_evaluate_function(DirectX::XMVECTOR x, float (*const fcn)(DirectX::XMFLOAT3 const &x, void *user_data), void *user_data, float const dirn);

static inline void internal_brx_fminsearch_3d_sort(DirectX::XMVECTOR inout_X[4], DirectX::XMVECTOR &inout_f);

static inline void brx_fminsearch_3d(float (*fcn)(DirectX::XMFLOAT3 const &x, void *user_data), void *user_data, DirectX::XMFLOAT3 &inout_x)
{
    //
    // Copyright (C) 2003-2025 The Octave Project Developers
    //
    // See the file COPYRIGHT.md in the top-level directory of this
    // distribution or <https://octave.org/copyright/>.
    //
    // This file is part of Octave.
    //
    // Octave is free software: you can redistribute it and/or modify it
    // under the terms of the GNU General Public License as published by
    // the Free Software Foundation, either version 3 of the License, or
    // (at your option) any later version.
    //
    // Octave is distributed in the hope that it will be useful, but
    // WITHOUT ANY WARRANTY; without even the implied warranty of
    // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    // GNU General Public License for more details.
    //
    // You should have received a copy of the GNU General Public License
    // along with Octave; see the file COPYING.  If not, see
    // <https://www.gnu.org/licenses/>.
    //

    // https://github.com/gnu-octave/octave/blob/default/scripts/optimization/fminsearch.m

    // References:
    // N. J. Higham, Optimization by direct search in matrix computations, SIAM J. Matrix Anal. Appl, 14(2): 317-333, 1993.
    // C. T. Kelley, Iterative Methods for Optimization, Society for Industrial and Applied Mathematics, Philadelphia, PA, 1999.

    constexpr uint32_t const n = 3U;
    constexpr uint32_t const max_iter = 200U * n;

    // Tolerance for cgce test based on relative size of simplex
    constexpr float const tol_x = 1E-4F;

    // Tolerance for cgce test based on step in function value.
    constexpr float const tol_fun = 1E-4F;

    //  Use function to minimize, not maximize
    constexpr float const dirn = -1.0F;

    DirectX::XMVECTOR X[n + 1U];
    DirectX::XMVECTOR f;
    {
        DirectX::XMVECTOR x0 = DirectX::XMLoadFloat3(&inout_x);
        float f0 = internal_brx_fminsearch_3d_guarded_evaluate_function(x0, fcn, user_data, dirn);

        // Set up initial simplex.
        float scale;
        {
            DirectX::XMVECTOR abs_x0_012 = DirectX::XMVectorAbs(x0);

            DirectX::XMFLOAT3 float3_abs_x0_012;
            DirectX::XMStoreFloat3(&float3_abs_x0_012, abs_x0_012);

            scale = std::max(1.0F, std::max(std::max(float3_abs_x0_012.x, float3_abs_x0_012.y), float3_abs_x0_012.z));
        }

        // Regular simplex - all edges have same length.
        // Generated from construction given in reference [18, pp. 80-81] of [1].
        DirectX::XMFLOAT2 alpha;
        DirectX::XMStoreFloat2(&alpha, DirectX::XMVectorScale(DirectX::XMVectorSet(std::sqrt(static_cast<float>(n) + 1.0F) - 1.0F + static_cast<float>(n), std::sqrt(static_cast<float>(n) + 1.0F) - 1.0F, 0.0F, 0.0F), scale / (static_cast<float>(n) * std::sqrt(2.0F))));

        static_assert(3U == n, "");
        DirectX::XMVECTOR x1 = DirectX::XMVectorAdd(x0, DirectX::XMVectorSet(alpha.x, alpha.y, alpha.y, 0.0F));
        float f1 = internal_brx_fminsearch_3d_guarded_evaluate_function(x1, fcn, user_data, dirn);
        DirectX::XMVECTOR x2 = DirectX::XMVectorAdd(x0, DirectX::XMVectorSet(alpha.y, alpha.x, alpha.y, 0.0F));
        float f2 = internal_brx_fminsearch_3d_guarded_evaluate_function(x2, fcn, user_data, dirn);
        DirectX::XMVECTOR x3 = DirectX::XMVectorAdd(x0, DirectX::XMVectorSet(alpha.y, alpha.y, alpha.x, 0.0F));
        float f3 = internal_brx_fminsearch_3d_guarded_evaluate_function(x3, fcn, user_data, dirn);

        X[0] = x0;
        X[1] = x1;
        X[2] = x2;
        X[3] = x3;

        f = DirectX::XMVectorSet(f0, f1, f2, f3);
    }

    // sort descend
    internal_brx_fminsearch_3d_sort(X, f);

    // reflect
    constexpr float const alpha = 1.0F;

    // expand
    constexpr float const gamma = 2.0F;

    // contract
    constexpr float const beta = 0.5F;

    // Outer (and only) loop.
    for (uint32_t k = 1U; k <= max_iter; ++k)
    {
        // assert(k < max_iter);

        // Stopping Test 3 - converged? The first part is test (4.3) in [1].
        float size_simplex;
        {
            DirectX::XMVECTOR splat_one = DirectX::XMVectorSplatOne();

            float norm_delta_x1 = DirectX::XMVectorGetX(DirectX::XMVector3Dot(DirectX::XMVectorAbs(DirectX::XMVectorSubtract(X[1], X[0])), splat_one));
            float norm_delta_x2 = DirectX::XMVectorGetX(DirectX::XMVector3Dot(DirectX::XMVectorAbs(DirectX::XMVectorSubtract(X[2], X[0])), splat_one));
            float norm_delta_x3 = DirectX::XMVectorGetX(DirectX::XMVector3Dot(DirectX::XMVectorAbs(DirectX::XMVectorSubtract(X[3], X[0])), splat_one));

            float norm_delta_x123 = std::max(std::max(norm_delta_x1, norm_delta_x2), norm_delta_x3);

            float norm_x0 = DirectX::XMVectorGetX(DirectX::XMVector3Dot(DirectX::XMVectorAbs(X[0]), splat_one));

            size_simplex = norm_delta_x123 / std::max(1.0F, norm_x0);
        }

        float step_f;
        {
            DirectX::XMVECTOR abs_delta_f0123 = DirectX::XMVectorAbs(DirectX::XMVectorSubtract(DirectX::XMVectorSplatX(f), f));

            DirectX::XMFLOAT4 float3_abs_delta_f0123;
            DirectX::XMStoreFloat4(&float3_abs_delta_f0123, abs_delta_f0123);

            step_f = std::max(std::max(float3_abs_delta_f0123.y, float3_abs_delta_f0123.z), float3_abs_delta_f0123.w);
        }

        if (size_simplex <= tol_x && step_f <= tol_fun)
        {
            // std::printf("Algorithm converged. Simplex size %9.4e <= %9.4e and step in function value %9.4e <= %9.4e\n", size_simplex, tol_x, step_f, tol_fun);
            break;
        }

        // One step of the Nelder-Mead simplex algorithm
        // NJH: Altered function calls and changed CNT to NF.
        // Changed each 'fr < f(1)' type test to '>' for maximization and re-ordered function values after sort.

        // Mean value
        static_assert(3U == n, "");
        DirectX::XMVECTOR xbar = DirectX::XMVectorScale(DirectX::XMVectorAdd(DirectX::XMVectorAdd(X[0], X[1]), X[2]), 1.0F / static_cast<float>(n));

        DirectX::XMVECTOR xr = DirectX::XMVectorLerp(xbar, X[n], -alpha);

        float fr = internal_brx_fminsearch_3d_guarded_evaluate_function(xr, fcn, user_data, dirn);

        // reflect
        DirectX::XMVECTOR xk = xr;
        float fk = fr;

        static_assert(3U == n, "");
        float f0;
        float f_n_minus_1;
        float f_n;
        {
            DirectX::XMFLOAT4 float4_f;
            DirectX::XMStoreFloat4(&float4_f, f);

            f0 = float4_f.x;
            f_n_minus_1 = float4_f.z;
            f_n = float4_f.w;
        }

        if (fr > f_n_minus_1)
        {
            if (fr > f0)
            {
                DirectX::XMVECTOR xe = DirectX::XMVectorLerp(xbar, xr, gamma);

                float fe = internal_brx_fminsearch_3d_guarded_evaluate_function(xe, fcn, user_data, dirn);

                if (fe > f0)
                {
                    // expand
                    xk = xe;
                    fk = fe;
                }
            }
        }
        else
        {
            DirectX::XMVECTOR xt = X[n];

            static_assert(3U == n, "");
            float ft = f_n;

            if (fr > ft)
            {
                xt = xr;
                ft = fr;
            }

            DirectX::XMVECTOR xc = DirectX::XMVectorLerp(xbar, xt, beta);

            float fc = internal_brx_fminsearch_3d_guarded_evaluate_function(xc, fcn, user_data, dirn);

            if (fc > f_n_minus_1)
            {
                // contract
                xk = xc;
                fk = fc;
            }
            else
            {

                static_assert(3U == n, "");
                DirectX::XMVECTOR x1 = DirectX::XMVectorScale(DirectX::XMVectorAdd(X[0], X[1]), 0.5F);
                float f1 = internal_brx_fminsearch_3d_guarded_evaluate_function(x1, fcn, user_data, dirn);
                DirectX::XMVECTOR x2 = DirectX::XMVectorScale(DirectX::XMVectorAdd(X[0], X[2]), 0.5F);
                float f2 = internal_brx_fminsearch_3d_guarded_evaluate_function(x2, fcn, user_data, dirn);

                X[1] = x1;
                X[2] = x2;

                f = DirectX::XMVectorSetZ(DirectX::XMVectorSetY(f, f1), f2);

                // shrink
                xk = DirectX::XMVectorScale(DirectX::XMVectorAdd(X[0], X[n]), 0.5F);
                fk = internal_brx_fminsearch_3d_guarded_evaluate_function(xk, fcn, user_data, dirn);
            }
        }

        X[n] = xk;
        static_assert(3U == n, "");
        f = DirectX::XMVectorSetW(f, fk);

        internal_brx_fminsearch_3d_sort(X, f);
    }
    // End of outer (and only) loop.

    // Finished
    DirectX::XMStoreFloat3(&inout_x, X[0]);
}

static inline float internal_brx_fminsearch_3d_guarded_evaluate_function(DirectX::XMVECTOR x, float (*const fcn)(DirectX::XMFLOAT3 const &x, void *user_data), void *user_data, float const dirn)
{
    DirectX::XMFLOAT3 float3_x;
    DirectX::XMStoreFloat3(&float3_x, x);
    float f = dirn * fcn(float3_x, user_data);
    if (!std::isnan(f))
    {
        // Do Nothing
    }
    else
    {
        assert(false);
        f = dirn * static_cast<float>(INFINITY);
    }
    return f;
}

static inline void internal_brx_fminsearch_3d_sort(DirectX::XMVECTOR inout_X[4], DirectX::XMVECTOR &inout_f)
{
    uint32_t uint4_idx[4] = {0U, 1U, 2U, 3U};

    DirectX::XMVECTOR idx = DirectX::XMLoadInt4(uint4_idx);

    // Layer 1: pairs (0,1) and (2,3) in parallel
    {
        DirectX::XMVECTOR p1 = DirectX::XMVectorSwizzle<0, 0, 2, 2>(inout_f);
        DirectX::XMVECTOR p2 = DirectX::XMVectorSwizzle<1, 1, 3, 3>(inout_f);
        DirectX::XMVECTOR ip1 = DirectX::XMVectorSwizzle<0, 0, 2, 2>(idx);
        DirectX::XMVECTOR ip2 = DirectX::XMVectorSwizzle<1, 1, 3, 3>(idx);

        // [cmp01, cmp01, cmp23, cmp23]
        DirectX::XMVECTOR cmp = DirectX::XMVectorGreater(p2, p1);

        DirectX::XMVECTOR hi = DirectX::XMVectorSelect(p1, p2, cmp);
        DirectX::XMVECTOR lo = DirectX::XMVectorSelect(p2, p1, cmp);

        DirectX::XMVECTOR idx_hi = DirectX::XMVectorSelect(ip1, ip2, cmp);
        DirectX::XMVECTOR idx_lo = DirectX::XMVectorSelect(ip2, ip1, cmp);

        // [max01, min01, max23, min23]
        DirectX::XMVECTOR m1010 = DirectX::XMVectorSelectControl(1, 0, 1, 0);

        inout_f = DirectX::XMVectorSelect(lo, hi, m1010);

        idx = DirectX::XMVectorSelect(idx_lo, idx_hi, m1010);
    }

    // Layer 2: pairs (0,2) and (1,3) in parallel
    {
        DirectX::XMVECTOR p1 = DirectX::XMVectorSwizzle<0, 1, 0, 1>(inout_f);
        DirectX::XMVECTOR p2 = DirectX::XMVectorSwizzle<2, 3, 2, 3>(inout_f);
        DirectX::XMVECTOR ip1 = DirectX::XMVectorSwizzle<0, 1, 0, 1>(idx);
        DirectX::XMVECTOR ip2 = DirectX::XMVectorSwizzle<2, 3, 2, 3>(idx);

        // [cmp02, cmp13, cmp02, cmp13]
        DirectX::XMVECTOR cmp = DirectX::XMVectorGreater(p2, p1);

        DirectX::XMVECTOR hi = DirectX::XMVectorSelect(p1, p2, cmp);
        DirectX::XMVECTOR lo = DirectX::XMVectorSelect(p2, p1, cmp);

        DirectX::XMVECTOR idx_hi = DirectX::XMVectorSelect(ip1, ip2, cmp);
        DirectX::XMVECTOR idx_lo = DirectX::XMVectorSelect(ip2, ip1, cmp);

        // [max02, max13, min02, min13]
        DirectX::XMVECTOR m1100 = DirectX::XMVectorSelectControl(1, 1, 0, 0);

        inout_f = DirectX::XMVectorSelect(lo, hi, m1100);

        idx = DirectX::XMVectorSelect(idx_lo, idx_hi, m1100);
    }

    // Layer 3: single pair (1,2)
    // pair (0,3) is guaranteed but comparing them in parallel will not break the correctness
    {
        DirectX::XMVECTOR p1 = DirectX::XMVectorSwizzle<0, 1, 1, 0>(inout_f);
        DirectX::XMVECTOR p2 = DirectX::XMVectorSwizzle<3, 2, 2, 3>(inout_f);
        DirectX::XMVECTOR ip1 = DirectX::XMVectorSwizzle<0, 1, 1, 0>(idx);
        DirectX::XMVECTOR ip2 = DirectX::XMVectorSwizzle<3, 2, 2, 3>(idx);

        // [cmp03, cmp12, cmp12, cmp03]
        DirectX::XMVECTOR cmp = DirectX::XMVectorGreater(p2, p1);

        DirectX::XMVECTOR hi = DirectX::XMVectorSelect(p1, p2, cmp);
        DirectX::XMVECTOR lo = DirectX::XMVectorSelect(p2, p1, cmp);

        DirectX::XMVECTOR idx_hi = DirectX::XMVectorSelect(ip1, ip2, cmp);
        DirectX::XMVECTOR idx_lo = DirectX::XMVectorSelect(ip2, ip1, cmp);

        // [max03, max12, min12, min03]
        DirectX::XMVECTOR m1100 = DirectX::XMVectorSelectControl(1, 1, 0, 0);

        inout_f = DirectX::XMVectorSelect(lo, hi, m1100);

        idx = DirectX::XMVectorSelect(idx_lo, idx_hi, m1100);
    }

    DirectX::XMStoreInt4(uint4_idx, idx);

    assert(uint4_idx[0] != uint4_idx[1]);
    assert(uint4_idx[0] != uint4_idx[2]);
    assert(uint4_idx[0] != uint4_idx[3]);
    assert(uint4_idx[1] != uint4_idx[2]);
    assert(uint4_idx[1] != uint4_idx[3]);
    assert(uint4_idx[2] != uint4_idx[3]);

    DirectX::XMVECTOR tmp_X[4];
    tmp_X[0] = inout_X[uint4_idx[0]];
    tmp_X[1] = inout_X[uint4_idx[1]];
    tmp_X[2] = inout_X[uint4_idx[2]];
    tmp_X[3] = inout_X[uint4_idx[3]];

    inout_X[0] = tmp_X[0];
    inout_X[1] = tmp_X[1];
    inout_X[2] = tmp_X[2];
    inout_X[3] = tmp_X[3];
}

#endif
