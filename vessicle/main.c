//
//  main.c
//  vessl
//
//  Created by Damien Di Fede on 10/20/15.
//  Copyright © 2015 compartmental. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vessl.h"

static int render_saw( float* buffer, unsigned int nBufferFrames, void* userData )
{
    float* t = (float*)userData;
    
    for( int i = 0; i < nBufferFrames; ++i )
    {
        buffer[i*2]   += *t * 0.01f;
        buffer[i*2+1] += *t * 0.01f;
    }
    
    *t = *t + 0.2f;
    
    if ( *t >= 1 )
    {
        *t = *t - 2;
    }
    
    return 0;
}

static int render_noise( float* buffer, unsigned int nBufferFrames, void* userData )
{
    float amp = *(float*)userData;
    
    for( int i = 0; i < nBufferFrames; ++i )
    {
        buffer[i*2]   += (float)rand()/RAND_MAX * amp;
        buffer[i*2+1] += (float)rand()/RAND_MAX * amp;
    }
    
    return 0;
}

static int render_mod( float* buffer, unsigned int nBufferFrames, void* userData )
{
    float mod = (float)rand()/RAND_MAX;
    
    for( int i = 0; i < nBufferFrames; ++i )
    {
        buffer[i*2]   *= mod;
        buffer[i*2+1] *= mod;
    }
    
    return 0;
}

static int render_unit_noise( float* sampleFrame, void* userData )
{
    float amp = *(float*)userData;
    
    sampleFrame[0] += (float)rand()/RAND_MAX * amp;
    sampleFrame[1] += (float)rand()/RAND_MAX * amp;
    
    return 0;
}

static float readBuffer[1024];

static int render_file( float* buffer, unsigned int nBufferFrames, void* userData )
{
    const char* filePath = (const char*)userData;
    
    vessl.e( "read", "file", filePath, &readBuffer, &nBufferFrames );
    
    for(unsigned int i = 0; i < nBufferFrames; ++i)
    {
        buffer[i*2] += readBuffer[i*2];
        buffer[i*2+1] += readBuffer[i*2+1];
    }
    
    return 0;
}

int main(int argc, const char * argv[])
{
    vessl.e( "hello" );
    
    vessl.e( "open", "output" );
    
    long bs = vessl.g( "output", "bs" ).l;
    printf("output has buffer size of %ld\n", bs);
    
    unsigned int sr = vessl.g( "output", "sr" ).u;
    printf("output has sample rate of %u\n", sr);
    
    /**
    float t = 0;
    vessl.e( "render", "output", "callback", &render_saw, &t );
    
    float amp = 0.01f;
    vessl.e( "render", "output", "callback", &render_noise, &amp );
    **/
    
    //float amp = 0.1f;
    //vessl.e( "render", "output", "unit", &render_unit_noise, &amp );
    
    const char * testFile = "/Users/ddf/dev/vessl/data/011513.wav";
    
    //vessl.e( "open", "file", testFile );
    
    //vessl.e( "render", "output", "callback", &render_file, testFile );
    
    //vessl.e( "play", testFile );
    
    //vessl.e( "render", "output", "callback", &render_mod );
    
    //vessl.e( "render", "output", "expression", "( t*8 * ((t*8>>12|t*8>>4) & 63 & t*8>>16)) * 32" );
    //vessl.e( "render", "output", "expression", "(((p/4) | (t*512)) & ((p/4) | (t*128)))" );
    
    // descending tone loop, take out the 50 - and it will be ascending
    //vessl.e( "render", "output", "expression", "t*(200 + (50 - m/50)%100)" );
    
    // a square wave with a fade in
    //vessl.e( "render", "output", "expression", "t*(m%2)" );
    
    // square wave with a weird polyrhythm
    //vessl.e( "render", "output", "expression", "t*(6 - m%2)" );
    
    // a really fucked up variation on the fade in square
    //vessl.e( "render", "output", "expression", "t*(m%((m*3/100)%12+2))<<2" );
    
    //vessl.e( "render", "output", "expression", "t*10 * ((t/441)+1)" );
    
    // this one gets really cool if you let it go for a while
    //vessl.e( "play", "expr", "(t<<t/(1024*8) | t>>t/16 & t>>t/32) / (t%(t/512+1) + 1) * 32" );
    
    // started as an attempt to control amplitude of a saw, got out of hand
    //vessl.e( "render", "output", "expression", "(r/2 - (256*(m/16%16)) + (t*(m/16%16)%(512*(m/16%16)+1))) * (m/16)" );
    
    // pure sine tone!
    //vessl.e( "play", "expr", "$(t*128)" );
    
    // frequency modulation of a saw using the sin operator.
    // the number inside the parens controls speed of modulation
    // the number following the parens controls frequency range modulated over
    //vessl.e( "play", "expr", "t*128 + $(t*5)*2" );
    
    // you can sort of amplitude modulate like this, but there is some extra stuff that gets in there
    // this does not sound like a saw mixed with a sine
    //vessl.e( "render", "output", "expression", "t*128 | $(t*5)" );
    
    // this is sort of surprising, shifting p by 10 gives a slightly glitchier version
    //vessl.e( "play", "expr", "t*128 | $(p>>12)*16" );
    
    // this is a cute little ditty
    //vessl.e( "play", "expr", "(t*128 + $(t)) | t>>(t%(8*r))/r | t>>128" );
    
    // this is a bit more aggressive, but has lots of nice texture
    //vessl.e( "render", "output", "expression", "(t*64 + $(t^$(m/2000))*$(m/2000)) | t*32" );
    
    //vessl.e( "render", "output", "expression", "(1 + $(m)%32) ^ (t*128 & t*64 & t*32) | (p/16)<<p%4 | $(p/128)>>p%4" );
    
    // this is strange but also not just noise or silence
    //vessl.e( "render", "output", "expression", "$t<<p^m-2*t|r+r|t*2-m^p>>$t" );
    
    // i dunno a weird sine based thing. kinda trying to write an expression similar to an ellipse equation
    //vessl.e( "render", "output", "expression", "(m/250+1)*$(t*128) | (m/500+1)*$((t+r/2*128))" );
    
    // based on equation for a hyperbola
    //vessl.e( "render", "output", "expression", "(m/128 * (t + 12345) + 32) | -(m/128 * (t + 12345) + 32)" );
    
    // a parabola?
    //vessl.e( "render", "output", "expression", "m*t*t | m*t | m" );
    
    // cumulative moving average, compare to the result of just using the x[n+1] part of it
    // this starts out pretty active and gradually becomes more stable, but doesn't ever seem to get stuck
    //vessl.e( "play", "expr", "p + ( ((t+1)*256 ^ (t+1)*64 & (t+1)*32) - p)/(t+1)");
    //vessl.e( "render", "output", "expression", "(t+1)*256 ^ (t+1)*64");
    
    // exponential moving average? not as interesting
    //vessl.e( "play", "expr", "$75 * (t*8>>12 | t*8>>8) | $(100-75) * p");
    
    vessl_value value;
    
    char* line = NULL;
    size_t linecap = 0;
    ssize_t linelen;
    while((linelen = getline(&line, &linecap, stdin)) > 0)
    {
        // parse the line into arguments for vessl, we will note what
        const char* arg1 = line;
        const char* arg2 = 0;
        const char* arg3 = 0;
        const char* arg4 = 0;
        for(int i = 0; i < linelen; ++i)
        {
            if ( !arg2 && line[i] == ' ' )
            {
                line[i] = '\0';
                arg2 = line + (i+1);
            }
            else if ( !arg3 && line[i] == ' ' )
            {
                line[i] = '\0';
                arg3 = line + (i+1);
            }
            else if ( !arg4 && line[i] == ' ' )
            {
                line[i] = '\0';
                arg4 = line + (i+1);
            }
        }
        
        // nuke the newline at the end of the string
        line[linelen-1] = '\0';
        
        size_t argLen = strlen(arg1);
        
        // exit the program
        if ( strncmp(arg1, "exit", argLen)==0 )
        {
            break;
        }
        // get a value
        else if ( strncmp(arg1, "get", argLen)==0 )
        {
            if ( arg2 == 0 || arg3 == 0 )
            {
                printf("[vessl] get syntax is: get objectName propertyName\n");
            }
            else
            {
                value = vessl.g( arg2, arg3 );
                printf( "%s = ", arg3 );
                switch (value.t)
                {
                    case VT_NONE:
                        printf("NULL\n");
                        break;
                    
                    case VT_DOUBLE:
                        printf("%f\n", value.d);
                        break;
                        
                    case VT_FLOAT:
                        printf("%f\n", value.f);
                        break;
                        
                    case VT_INT:
                        printf("%d\n", value.i);
                        break;
                        
                    case VT_LONG:
                        printf("%ld\n", value.l);
                        break;
                        
                    case VT_UINT:
                        printf("%u\n", value.u);
                        break;
                }
            }
        }
        else if ( strncmp(arg1, "set", argLen)==0 )
        {
            if ( arg2 == 0 || arg3 == 0 || arg4 == 0 )
            {
                printf("set syntax is: set objectName propertyName value\n");
            }
            else
            {
                char* end;
                long num = strtol(arg4, &end, 10);
                if ( end != arg4 )
                {
                    value.u = (unsigned int)num;
                    value.t = VT_UINT;
                    vessl.s( arg2, arg3, value );
                }
                else
                {
                    printf("[vessl] could not set %s because converting %s to a number did not work\n", arg3, arg4);
                }
            }
        }
        // do a define
        else if ( strncmp(arg1, "define", argLen)==0 )
        {
            vessl.d( arg2, arg3 );
        }
        else
        {
            vessl.e( arg1, arg2, arg3 );
        }
    }
    
    vessl.e( "close", "file", testFile );
    
    vessl.e( "close", "output" );
    
    return 0;
}
