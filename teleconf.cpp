/*
 * File: teleconf.cpp
 */

#include "usermodfx.h"
#include "biquad.hpp"

// filter parameters 
#define LPF_FC 1400.f
#define HPF_FC 300.f
#define LPF_Q  0.707f

// filters before downsampling
static dsp::BiQuad s_bq_lpf;
static dsp::BiQuad s_bqs_lpf;
static dsp::BiQuad s_bq_hpf;
static dsp::BiQuad s_bqs_hpf;
// filters after downsampling
static dsp::BiQuad s_bq_lpf_out;
static dsp::BiQuad s_bqs_lpf_out;
static dsp::BiQuad s_bq_lpf_out2;
static dsp::BiQuad s_bqs_lpf_out2;

static const float s_fs_recip = 1.f / 48000.f;

static float dry = 0.f;
static float wet = 1.f;
static uint32_t count = 0;
static float lastmy;
static float lastsy;

void init_lpf(const float f, const float q) {
    float wc = s_bq_lpf.mCoeffs.wc(f, s_fs_recip);
    s_bq_lpf.mCoeffs.setSOLP(fx_tanpif(wc), q);
    s_bqs_lpf.mCoeffs     = s_bq_lpf.mCoeffs;
    s_bq_lpf_out.mCoeffs  = s_bq_lpf.mCoeffs;
    s_bqs_lpf_out.mCoeffs = s_bq_lpf.mCoeffs;
    s_bq_lpf_out2.mCoeffs  = s_bq_lpf.mCoeffs;
    s_bqs_lpf_out2.mCoeffs = s_bq_lpf.mCoeffs;
}

void MODFX_INIT(uint32_t platform, uint32_t api)
{
    s_bq_lpf.flush();
    s_bqs_lpf.flush();
    s_bq_lpf_out.flush();
    s_bqs_lpf_out.flush();
    s_bq_lpf_out2.flush();
    s_bqs_lpf_out2.flush();
    init_lpf(LPF_FC, LPF_Q);
    float wc = s_bq_hpf.mCoeffs.wc(HPF_FC, s_fs_recip);
    s_bq_hpf.flush();
    s_bqs_hpf.flush();
    s_bq_hpf.mCoeffs.setFOHP(fx_tanpif(wc));
    s_bqs_hpf.mCoeffs = s_bq_hpf.mCoeffs;
}

__fast_inline float g711(const float s) {
    q15_t val = f32_to_q15(s);
    int16_t sign = (val < 0) ? -1 : 1;
    val = q15abs(val);

    uint16_t mask = 1 << 14;
    int i;
    for(i = 0; i < 6; i++) {
        if (val & mask)
            break;
        else
            val <<=1;
    }
    val &= 0x7c00;
    val >>= i;
    val = val * sign;
    return q15_to_f32(val);
}

void MODFX_PROCESS(const float *main_xn, float *main_yn,
                   const float *sub_xn,  float *sub_yn,
                   uint32_t frames)
{
  const float * mx = main_xn;
  float * __restrict my = main_yn;
  const float * my_e = my + 2*frames;

  const float * sx = sub_xn;
  float * __restrict sy = sub_yn;
  
  float vmx;
  float vsx;
  for (; my != my_e; ) {
      vmx = s_bq_hpf.process_fo(s_bq_lpf.process_so(*mx));
      vsx = s_bqs_hpf.process_fo(s_bqs_lpf.process_so(*sx));

      if (count == 0) {
          lastmy = g711(vmx);
          lastsy = g711(vsx);
      }
      count = (count + 1) % 6;

      *my++ = dry * (*mx++) + wet * \
          s_bq_lpf_out2.process_so(s_bq_lpf_out.process_so(lastmy));
      *sy++ = dry * (*sx++) + wet * \
          s_bq_lpf_out2.process_so(s_bqs_lpf_out.process_so(lastsy));
  }
}

void MODFX_PARAM(uint8_t index, int32_t value)
{
  const float valf = q31_to_f32(value);

  switch (index) {
  case k_user_modfx_param_time:
      init_lpf(LPF_FC, LPF_Q  + 1.6 * valf);
    break;
  case k_user_modfx_param_depth:
      dry = valf;
      wet = 1.0 - dry;
    break;
  default:
    break;
  }
}
