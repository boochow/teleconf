/*
 * File: teleconf.cpp
 */

#include "usermodfx.h"
#include "biquad.hpp"

// filter parameters 
#define LPF_FC 3000.f
#define HPF_FC 300.f
#define LPF_Q  0.707f

// filters before downsampling
static dsp::BiQuad s_bq_lpf_r;
static dsp::BiQuad s_bqs_lpf_r;
static dsp::BiQuad s_bq_hpf_r;
static dsp::BiQuad s_bqs_hpf_r;
static dsp::BiQuad s_bq_lpf_l;
static dsp::BiQuad s_bqs_lpf_l;
static dsp::BiQuad s_bq_hpf_l;
static dsp::BiQuad s_bqs_hpf_l;
// filters after downsampling
static dsp::BiQuad s_bq_lpf_out_r;
static dsp::BiQuad s_bqs_lpf_out_r;
static dsp::BiQuad s_bq_lpf_out2_r;
static dsp::BiQuad s_bqs_lpf_out2_r;
static dsp::BiQuad s_bq_lpf_out_l;
static dsp::BiQuad s_bqs_lpf_out_l;
static dsp::BiQuad s_bq_lpf_out2_l;
static dsp::BiQuad s_bqs_lpf_out2_l;

static const float s_fs_recip = 1.f / 48000.f;

static float dry = 0.f;
static float wet = 1.f;
static uint32_t count = 0;
static float lastmy_r;
static float lastsy_r;
static float lastmy_l;
static float lastsy_l;

void init_lpf(const float f, const float q) {
    float wc = s_bq_lpf_r.mCoeffs.wc(f, s_fs_recip);
    s_bq_lpf_r.mCoeffs.setSOLP(fx_tanpif(wc), q);
    s_bqs_lpf_r.mCoeffs     = s_bq_lpf_r.mCoeffs;
    s_bq_lpf_out_r.mCoeffs  = s_bq_lpf_r.mCoeffs;
    s_bqs_lpf_out_r.mCoeffs = s_bq_lpf_r.mCoeffs;
    s_bq_lpf_out2_r.mCoeffs  = s_bq_lpf_r.mCoeffs;
    s_bqs_lpf_out2_r.mCoeffs = s_bq_lpf_r.mCoeffs;

    s_bq_lpf_l.mCoeffs     = s_bq_lpf_r.mCoeffs;
    s_bqs_lpf_l.mCoeffs     = s_bq_lpf_r.mCoeffs;
    s_bq_lpf_out_l.mCoeffs  = s_bq_lpf_r.mCoeffs;
    s_bqs_lpf_out_l.mCoeffs = s_bq_lpf_r.mCoeffs;
    s_bq_lpf_out2_l.mCoeffs  = s_bq_lpf_r.mCoeffs;
    s_bqs_lpf_out2_l.mCoeffs = s_bq_lpf_r.mCoeffs;
}

void MODFX_INIT(uint32_t platform, uint32_t api)
{
    s_bq_lpf_r.flush();
    s_bqs_lpf_r.flush();
    s_bq_lpf_out_r.flush();
    s_bqs_lpf_out_r.flush();
    s_bq_lpf_out2_r.flush();
    s_bqs_lpf_out2_r.flush();

    s_bq_lpf_l.flush();
    s_bqs_lpf_l.flush();
    s_bq_lpf_out_l.flush();
    s_bqs_lpf_out_l.flush();
    s_bq_lpf_out2_l.flush();
    s_bqs_lpf_out2_l.flush();

    init_lpf(LPF_FC, LPF_Q);

    float wc = s_bq_hpf_r.mCoeffs.wc(HPF_FC, s_fs_recip);
    s_bq_hpf_r.flush();
    s_bqs_hpf_r.flush();
    s_bq_hpf_l.flush();
    s_bqs_hpf_l.flush();

    s_bq_hpf_r.mCoeffs.setSOHP(fx_tanpif(wc), 0.5);
    s_bqs_hpf_r.mCoeffs = s_bq_hpf_r.mCoeffs;
    s_bq_hpf_l.mCoeffs = s_bq_hpf_r.mCoeffs;
    s_bqs_hpf_l.mCoeffs = s_bq_hpf_r.mCoeffs;
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

      // Left channel
      vmx = s_bq_hpf_l.process_so(s_bq_lpf_l.process_so(*mx));
      vsx = s_bqs_hpf_l.process_so(s_bqs_lpf_l.process_so(*sx));

      if (count == 0) {
          lastmy_l = g711(vmx);
          lastsy_l = g711(vsx);
      }

      *my++ = dry * (*mx++) + wet * \
          s_bq_lpf_out2_l.process_so(s_bq_lpf_out_l.process_so(lastmy_l));
      *sy++ = dry * (*sx++) + wet * \
          s_bq_lpf_out2_l.process_so(s_bqs_lpf_out_l.process_so(lastsy_l));

      // Right channel
      vmx = s_bq_hpf_r.process_so(s_bq_lpf_r.process_so(*mx));
      vsx = s_bqs_hpf_r.process_so(s_bqs_lpf_r.process_so(*sx));

      if (count == 0) {
          lastmy_r = g711(vmx);
          lastsy_r = g711(vsx);
      }

      *my++ = dry * (*mx++) + wet * \
          s_bq_lpf_out2_r.process_so(s_bq_lpf_out_r.process_so(lastmy_r));
      *sy++ = dry * (*sx++) + wet * \
          s_bq_lpf_out2_r.process_so(s_bqs_lpf_out_r.process_so(lastsy_r));


      count = (count + 1) % 6;
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
      wet = valf;
      dry = 1.0 - wet;
    break;
  default:
    break;
  }
}
