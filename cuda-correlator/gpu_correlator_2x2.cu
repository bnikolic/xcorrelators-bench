/*
Copyright (C) 2009 Rob van Nieuwpoort & John Romein
Astron
P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, nieuwpoort@astron.nl

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
__global__ void correlate_2x2(float4* devSamples, float *devVisibilities, unsigned nrTimes, unsigned nrTimesWidth, 
			      unsigned nrStations, unsigned nrChannels, unsigned nrCells,
			      unsigned nrBlocks, unsigned nrThreads, unsigned loopCount)
{
    unsigned myBlock   = blockIdx.x + gridDim.x * blockIdx.y;
    unsigned myThread  = threadIdx.x + blockDim.x * threadIdx.y;

    for (unsigned channel = myBlock; channel < nrChannels; channel += nrBlocks) {
	PREFETCH();
	for (unsigned cell = myThread; cell < loopCount; cell += nrThreads) {

	    unsigned stat0 = cellToStatX[cell];
	    unsigned stat2 = cellToStatY[cell];

#if !USE_ALTERNATE_MEMORY_LAYOUT
	    unsigned index0 = SAMPLE_INDEX(stat0, channel, 0, 0) / 2;
	    unsigned index2 = SAMPLE_INDEX(stat2, channel, 0, 0) / 2;
#else
	    unsigned index = SAMPLE_INDEX(0, channel, 0, 0) / 2;
#endif

	    float v0x2xr = 0, v0x2xi = 0;
	    float v0x2yr = 0, v0x2yi = 0;
	    float v0y2xr = 0, v0y2xi = 0;
	    float v0y2yr = 0, v0y2yi = 0;
	    float v1x2xr = 0, v1x2xi = 0;
	    float v1x2yr = 0, v1x2yi = 0;
	    float v1y2xr = 0, v1y2xi = 0;
	    float v1y2yr = 0, v1y2yi = 0;
	    float v0x3xr = 0, v0x3xi = 0;
	    float v0x3yr = 0, v0x3yi = 0;
	    float v0y3xr = 0, v0y3xi = 0;
	    float v0y3yr = 0, v0y3yi = 0;
	    float v1x3xr = 0, v1x3xi = 0;
	    float v1x3yr = 0, v1x3yi = 0;
	    float v1y3xr = 0, v1y3xi = 0;
	    float v1y3yr = 0, v1y3yi = 0;

	    for (unsigned time = 0; time < nrTimes; time ++) {
		float4 sample0 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample1 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample2 = {0.0f, 0.0f, 0.0f, 0.0f};
		float4 sample3 = {0.0f, 0.0f, 0.0f, 0.0f};

#if !DO_LOADS
		// A dummy condition that always evaluates to false.
		// Needed to trick the compiler
		if(nrTimes > 100000) {
#endif
		    // Just load the samples, even though cell >= nrCells.
		    // We are loading junk in that case, but still this is faster than testing.
#if USE_ALTERNATE_MEMORY_LAYOUT
#if USE_TEXTURE_CACHE
		    sample0 = tex1Dfetch(tex0, index + stat0+0);
		    sample1 = tex1Dfetch(tex0, index + stat0+1);

		    sample2 = tex1Dfetch(tex0, index + stat2+0);
		    sample3 = tex1Dfetch(tex0, index + stat2+1);
#else
		    sample0 = devSamples[index + stat0+0];
		    sample1 = devSamples[index + stat0+1];

		    sample2 = devSamples[index + stat2+0];
		    sample3 = devSamples[index + stat2+1];
#endif // USE_TEXTURE_CACHE
		    index += nrStations; // increase time
#else // !USE_ALTERNATE_MEMORY_LAYOUT
#if USE_TEXTURE_CACHE
		    sample0 = tex1Dfetch(tex0, index0 + time + 0 * nrTimesWidth);
		    sample1 = tex1Dfetch(tex0, index0 + time + 1 * nrTimesWidth);

		    sample2 = tex1Dfetch(tex0, index2 + time + 0 * nrTimesWidth);
		    sample3 = tex1Dfetch(tex0, index2 + time + 1 * nrTimesWidth);
#else
		    sample0 = devSamples[index0 + time + 0 * nrTimesWidth];
		    sample1 = devSamples[index0 + time + 1 * nrTimesWidth];

		    sample2 = devSamples[index2 + time + 0 * nrTimesWidth];
		    sample3 = devSamples[index2 + time + 1 * nrTimesWidth];
#endif // USE_TEXTURE_CACHE
#endif // USE_ALTERNATE_MEMORY_LAYOUT
#if !DO_LOADS
		}
#endif

#if SYNC_THREADS
		__syncthreads();
#endif

		v0x2xr  += sample0.x * sample2.x;
		v0x2xi  += sample0.y * sample2.x;
		v0x2yr  += sample0.x * sample2.z;
		v0x2yi  += sample0.y * sample2.z;
		v0y2xr  += sample0.z * sample2.x;
		v0y2xi  += sample0.w * sample2.x;
		v0y2yr  += sample0.z * sample2.z;
		v0y2yi  += sample0.w * sample2.z;

		v0x2xr  += sample0.y * sample2.y;
		v0x2xi  -= sample0.x * sample2.y;
		v0x2yr  += sample0.y * sample2.w;
		v0x2yi  -= sample0.x * sample2.w;
		v0y2xr  += sample0.w * sample2.y;
		v0y2xi  -= sample0.z * sample2.y;
		v0y2yr  += sample0.w * sample2.w;
		v0y2yi  -= sample0.z * sample2.w;

		// the load of sample 1 can be done here, but for
		// some reason it was slower...

		v1x2xr  += sample1.x * sample2.x;
		v1x2xi  += sample1.y * sample2.x;
		v1x2yr  += sample1.x * sample2.z;
		v1x2yi  += sample1.y * sample2.z;
		v1y2xr  += sample1.z * sample2.x;
		v1y2xi  += sample1.w * sample2.x;
		v1y2yr  += sample1.z * sample2.z;
		v1y2yi  += sample1.w * sample2.z;

		v1x2xr  += sample1.y * sample2.y;
		v1x2xi  -= sample1.x * sample2.y;
		v1x2yr  += sample1.y * sample2.w;
		v1x2yi  -= sample1.x * sample2.w;
		v1y2xr  += sample1.w * sample2.y;
		v1y2xi  -= sample1.z * sample2.y;
		v1y2yr  += sample1.w * sample2.w;
		v1y2yi  -= sample1.z * sample2.w;

		// load of sample3 can be done here, reusing sample2's registers

		v0x3xr  += sample0.x * sample3.x;
		v0x3xi  += sample0.y * sample3.x;
		v0x3yr  += sample0.x * sample3.z;
		v0x3yi  += sample0.y * sample3.z;
		v0y3xr  += sample0.z * sample3.x;
		v0y3xi  += sample0.w * sample3.x;
		v0y3yr  += sample0.z * sample3.z;
		v0y3yi  += sample0.w * sample3.z;

		v0x3xr  += sample0.y * sample3.y;
		v0x3xi  -= sample0.x * sample3.y;
		v0x3yr  += sample0.y * sample3.w;
		v0x3yi  -= sample0.x * sample3.w;
		v0y3xr  += sample0.w * sample3.y;
		v0y3xi  -= sample0.z * sample3.y;
		v0y3yr  += sample0.w * sample3.w;
		v0y3yi  -= sample0.z * sample3.w;


		v1x3xr  += sample1.x * sample3.x;
		v1x3xi  += sample1.y * sample3.x;
		v1x3yr  += sample1.x * sample3.z;
		v1x3yi  += sample1.y * sample3.z;
		v1y3xr  += sample1.z * sample3.x;
		v1y3xi  += sample1.w * sample3.x;
		v1y3yr  += sample1.z * sample3.z;
		v1y3yi  += sample1.w * sample3.z;

		v1x3xr  += sample1.y * sample3.y;
		v1x3xi  -= sample1.x * sample3.y;
		v1x3yr  += sample1.y * sample3.w;
		v1x3yi  -= sample1.x * sample3.w;
		v1y3xr  += sample1.w * sample3.y;
		v1y3xi  -= sample1.z * sample3.y;
		v1y3yr  += sample1.w * sample3.w;
		v1y3yi  -= sample1.z * sample3.w;
	    }

#if !DO_STORES
	    if(nrTimes > 100000) { // Dummy condition that always evaluates to false; trick the compiler
#endif
		if (cell < nrCells) {
		    unsigned baseline = BASELINE(stat0, stat2);
		    float4 *dst = (float4*)(devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		    dst[0].x = v0x2xr;
		    dst[0].y = v0x2xi;
		    dst[0].z = v0x2yr;
		    dst[0].w = v0x2yi;
		    dst[1].x = v0y2xr;
		    dst[1].y = v0y2xi;
		    dst[1].z = v0y2yr;
		    dst[1].w = v0y2yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v1x2xr;
		    dst[0].y = v1x2xi;
		    dst[0].z = v1x2yr;
		    dst[0].w = v1x2yi;
		    dst[1].x = v1y2xr;
		    dst[1].y = v1y2xi;
		    dst[1].z = v1y2yr;
		    dst[1].w = v1y2yi;

		    baseline = BASELINE(stat0, stat2+1);
		    dst = (float4*)(devVisibilities + 2 * VISIBILITIES_INDEX(baseline, channel, 0, 0));

		    dst[0].x = v0x3xr;
		    dst[0].y = v0x3xi;
		    dst[0].z = v0x3yr;
		    dst[0].w = v0x3yi;
		    dst[1].x = v0y3xr;
		    dst[1].y = v0y3xi;
		    dst[1].z = v0y3yr;
		    dst[1].w = v0y3yi;

		    // increase baseline by one
		    dst += nrChannels * 2;

		    dst[0].x = v1x3xr;
		    dst[0].y = v1x3xi;
		    dst[0].z = v1x3yr;
		    dst[0].w = v1x3yi;
		    dst[1].x = v1y3xr;
		    dst[1].y = v1y3xi;
		    dst[1].z = v1y3yr;
		    dst[1].w = v1y3yi;
		}
#if !DO_STORES
	    }
#endif
	}
    }
}
