//
//  vessl.c
//  vessl
//
//  Created by Damien Di Fede on 10/19/15.
//  Copyright © 2015 compartmental. All rights reserved.
//

#include <stdio.h>
#include <map>
#include <list>
#include <limits.h>
#include <sndfile.h>
#include <RtAudio.h>
#include "vessl.h"
#include "expr_eval.hpp"

// forward declare so we can use it
const char* resolve_define(const char*);

enum render_type
{
    eCallback,
    eUnit
};

typedef struct
{
    render_type           type;
    vessl_render_callback callback;
    vessl_render_unit     unit;
    void*                 userData;
} render_callback;

typedef std::list<render_callback> render_callback_list;

static RtAudio              vessl_out;
static const unsigned int   vessl_out_default_bufferFrames = 1024;
static unsigned int         vessl_out_bufferFrames = vessl_out_default_bufferFrames;
static const unsigned int   vessl_out_default_sampleRate = 44100;
static unsigned int         vessl_out_sampleRate = vessl_out_default_sampleRate;
static render_callback_list vessl_out_render_callbacks;

static int output_render_callback( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames, double streamTime, RtAudioStreamStatus status, void *userData )
{
    float* buffer = (float*)outputBuffer;
    
    memset(buffer, 0, sizeof(float)*nBufferFrames*2);
    
    for( auto iter = vessl_out_render_callbacks.begin(); iter != vessl_out_render_callbacks.end(); ++iter )
    {
        switch( iter->type )
        {
            case eCallback:
            {
                if ( iter->callback( buffer, nBufferFrames, iter->userData ) )
                {
                    iter = vessl_out_render_callbacks.erase(iter);
                }
            }
            break;
                
            case eUnit:
            {
                for(int i = 0; i < nBufferFrames; ++i)
                {
                    if ( iter->unit( buffer + (i*2), iter->userData ) )
                    {
                        iter = vessl_out_render_callbacks.erase(iter);
                        break;
                    }
                }
            }
            break;
        }
    }
    
    return 0;
}

static void open_output(va_list args)
{
    RtAudio::StreamParameters parameters;
    parameters.deviceId = vessl_out.getDefaultOutputDevice();
    parameters.nChannels = 2;
    parameters.firstChannel = 0;
    
    RtAudio::StreamOptions options;
    if ( vessl_out_bufferFrames == 0 )
    {
        options.flags |= RTAUDIO_MINIMIZE_LATENCY;
    }
    options.numberOfBuffers = 2;
    
    parameters.deviceId = vessl_out.getDefaultOutputDevice();
    try
    {
        vessl_out.openStream( &parameters, NULL, RTAUDIO_FLOAT32, vessl_out_sampleRate, &vessl_out_bufferFrames, &output_render_callback, (void*)0, &options );
    }
    catch( RtAudioError& e )
    {
        printf("[vessl] FAILED TO OPEN OUTPUT: %s\n", e.getMessage().c_str());
    }
    
    if ( vessl_out.isStreamOpen() )
    {
        try
        {
            vessl_out.startStream();
        }
        catch( RtAudioError& e )
        {
            printf("[vessl] FAILED TO START OUTPUT: %s\n", e.getMessage().c_str());
        }
    }
}

static void close_output(va_list args)
{
    try
    {
        // Stop the stream
        vessl_out.stopStream();
        
        if ( vessl_out.isStreamOpen() )
        {
            vessl_out.closeStream();
            
            printf("[vessl] output was closed.\n");
        }
    }
    catch (RtAudioError& e)
    {
        e.printMessage();
    }
}

static void render_output(va_list args)
{
    const char *          renderType = va_arg(args, const char *);
    if ( strcmp(renderType, "callback")==0 )
    {
        vessl_render_callback callback = va_arg(args, vessl_render_callback);
        void*                 userData = va_arg(args, void*);
    
        vessl_out_render_callbacks.push_back( { eCallback, callback, 0, userData } );
    }
    else if ( strcmp(renderType, "unit")==0 )
    {
        vessl_render_unit unit = va_arg(args, vessl_render_unit);
        void*             userData = va_arg(args, void*);
        
        vessl_out_render_callbacks.push_back( { eUnit, 0, unit, userData } );
    }
}

typedef enum
{
    ePlaying = 0,
    eOpen,
    eClosed
} vessl_file_state ;

typedef struct
{
    std::string         path;
    SNDFILE*            file;
    SF_INFO*            info;
    vessl_file_state    state;
} vessl_file;

static std::map<std::string, vessl_file*> vessl_files;

vessl_file* get_file(const char* named)
{
    named = resolve_define(named);
    auto iter = vessl_files.find(named);
    if ( iter != vessl_files.end() )
    {
        return iter->second;
    }
    
    return 0;
}

static void open_file(va_list args)
{
    const char * filePath = resolve_define( va_arg(args, const char*) );

    SF_INFO* info = new SF_INFO;
    info->format = 0;
    SNDFILE* file = sf_open(filePath, SFM_READ, info);
    if( file )
    {
        vessl_files[filePath] = new vessl_file { filePath, file, info, eOpen };
    }
    else
    {
        delete info;
        printf("[vessl] couldn't open %s because %s\n", filePath, sf_strerror(file));
    }
}

static void close_filep( vessl_file* handle )
{
    sf_close(handle->file);
    delete handle->info;
    delete handle;
}

static void close_file(va_list args)
{
    const char * filePath = va_arg(args, const char*);
    
    if ( vessl_file* file = get_file(filePath) )
    {
        if ( file->state == eOpen )
        {
            close_filep(file);
        }
        else
        {
            file->state = eClosed;
        }
        vessl_files.erase(vessl_files.find(file->path));
    }
}

static void read_file(va_list args)
{
    const char * filePath = va_arg(args, const char*);
    
    if ( vessl_file* vfile = get_file(filePath) )
    {
        float* buffer              = va_arg(args, float*);
        unsigned int* framesToread = va_arg(args, unsigned int*);
        sf_count_t count = *framesToread;
        *framesToread = (unsigned int)sf_readf_float(vfile->file, buffer, count);
    }
}

// reasonably large array to read from files
static float play_file_read_buffer[vessl_out_default_bufferFrames*4];

static int play_file_callback( float* buffer, unsigned int nBufferFrames, void* userData )
{
    vessl_file* handle = (vessl_file*)userData;
    int result = 0;
    
    switch( handle->state )
    {
        case ePlaying:
        {
            sf_count_t count = nBufferFrames;
            count = sf_readf_float(handle->file, play_file_read_buffer, count);
            
            switch( handle->info->channels )
            {
                case 1:
                    for(sf_count_t i = 0; i < count; ++i)
                    {
                        buffer[i*2]   += play_file_read_buffer[i];
                        buffer[i*2+1] += play_file_read_buffer[i];
                    }
                    break;
                    
                case 2:
                    count *= 2;
                    for(sf_count_t i = 0; i < count; ++i)
                    {
                        buffer[i] += play_file_read_buffer[i];
                    }
                    break;
            }
        }
        break;
            
        case eClosed:
        {
            close_filep(handle);
            result = 1;
        }
        break;
            
        case eOpen:
        {
            result = 1;
        }
        break;
    }
    
    return result;
}

static void play_file(va_list args)
{
    const char * filePath = resolve_define( va_arg(args, const char*) );
    
    if ( vessl_out.isStreamOpen() == false )
    {
        vessl.e( "open", "output" );
    }
    
    auto iter = vessl_files.find(filePath);
    
    // not open yet, so open it
    if ( iter == vessl_files.end() )
    {
        vessl.e( "open", "file", filePath );
        // make sure it opened properly
        iter = vessl_files.find(filePath);
    }

    // if it is open, start playing if it isn't playing already!
    if ( iter != vessl_files.end() && iter->second->state != ePlaying )
    {
        iter->second->state = ePlaying;
        vessl.e( "render", "output", "callback", &play_file_callback, iter->second );
    }
}

static void pause_file(va_list args)
{
    const char * filePath = va_arg(args, const char*);

    if ( vessl_file* file = get_file(filePath) )
    {
        file->state = eOpen;
    }
    else
    {
        printf("[vessl] could not pause %s because it is not open\nd", filePath);
    }
}

struct vessl_expression
{
    std::string  exp;
    float        amp;
    enum { PLAY, PAUSE, STOP } state;
    ExprEval     eval;
};

// the expression string is used as the key
static std::map<std::string, vessl_expression*> vessl_expressions;

static vessl_expression* get_expr(const char* named)
{
    named = resolve_define(named);
    auto iter = vessl_expressions.find(named);
    if ( iter == vessl_expressions.end() )
    {
        return 0;
    }
    
    return iter->second;
}

int vessl_expression_unit( float* sampleFrame, void* userData )
{
    vessl_expression* vexp = (vessl_expression*)userData;
    
    if ( vexp->state == vessl_expression::PLAY )
    {
        char* exp = const_cast<char*>(vexp->exp.c_str());
        
        unsigned int result = vexp->eval.Eval(exp);
        
        switch ( vexp->eval.GetErr() )
        {
            case EEE_DIVIDE_BY_ZERO:
            case EEE_PARENTHESIS:
            case EEE_WRONG_CHAR:
            {
                vexp->state = vessl_expression::STOP;
                printf("[vessl] error evaluating %s at %s: %d\n", exp, vexp->eval.GetErrPos(), vexp->eval.GetErr());
            }
                
            default: break;
        }
        
        if ( vexp->state == vessl_expression::STOP )
        {
            vessl_expressions.erase(vexp->exp);
            delete vexp;
            return 1;
        }
        
        vexp->eval.SetVar('p', result);
        
        unsigned int t = vexp->eval.GetVar('t') + 1;
        vexp->eval.SetVar('t', t);
        
        // milliseconds
        vexp->eval.SetVar('m', t/(vessl_out_sampleRate/1000));
        
        const unsigned int range = vexp->eval.GetVar('r');
        
        const float h = range / 2;
        
        const float sample = vexp->amp * (((float)(result%range) - h) / h);
        
        sampleFrame[0] += sample;
        sampleFrame[1] += sample;
    }
    
    return 0;
}

static vessl_expression* create_expr(const char* exp)
{
    exp = resolve_define(exp);
    
    vessl_expression* vexp = new vessl_expression { exp, 0.1f, vessl_expression::PLAY };
    // time, increased by 1 each time this is evaulated
    vexp->eval.SetVar('t', 0);
    // value of the previous evaluation
    vexp->eval.SetVar('p', 0);
    // milliseconds elapsed
    vexp->eval.SetVar('m', 0);
    // range that the eval output is wrapped to
    vexp->eval.SetVar('r', 1<<15);
    
    vessl_expressions[exp] = vexp;
    vessl_out_render_callbacks.push_back( { eUnit, 0, &vessl_expression_unit, vexp } );
    
    return vexp;
}

static void play_expr(va_list args)
{
    const char* exp = va_arg(args, const char*);
    if ( vessl_expression* vexp = get_expr(exp) )
    {
        vexp->state = vessl_expression::PLAY;
    }
    else
    {
        create_expr(exp);
    }
}

static void pause_expr(va_list args)
{
    const char* exp = va_arg(args, const char*);
    if ( vessl_expression* vexp = get_expr(exp) )
    {
        vexp->state = vessl_expression::PAUSE;
    }
}

//-- ACTIONS --------------------------------------------------

typedef void(* const vessl_action)(va_list);

struct StrCompare : public std::binary_function<const char*, const char*, bool>
{
public:
    bool operator() (const char* str1, const char* str2) const
    { return std::strcmp(str1, str2) < 0; }
};

typedef std::map<const char*, vessl_action, StrCompare> vessl_action_table;

static void call_action( vessl_action_table& table, va_list args )
{
    const char * name = va_arg(args, const char*);
    auto iter = table.find(name);
    if ( iter != table.end() )
    {
        iter->second(args);
    }
    else
    {
        printf("[vessl] i don't know what is %s\n", name);
    }
}

static vessl_action_table open_actions = {
    { "output", open_output },
    { "file", open_file }
};

static void open(va_list args)
{
    call_action( open_actions, args );
}

static vessl_action_table close_actions = {
    { "output", close_output },
    { "file", close_file }
};

static void close(va_list args)
{
    call_action( close_actions, args );
}

static vessl_action_table render_actions = {
    { "output", render_output }
};

static void render(va_list args)
{
    call_action( render_actions, args );
}

static vessl_action_table read_actions = {
    { "file", read_file }
};

static void read(va_list args)
{
    call_action( read_actions, args );
}

static vessl_action_table play_actions = {
    { "file", play_file },
    { "expr", play_expr }
};

static void play(va_list args)
{
    call_action(play_actions, args);
}

static vessl_action_table pause_actions = {
    { "file", pause_file },
    { "expr", pause_expr }
};

static void pause(va_list args)
{
    call_action(pause_actions, args);
}

static vessl_action_table vessl_actions = {
    { "open",   open },
    { "close",  close },
    { "render", render },
    { "read",   read },
    { "play",   play },
    { "pause",  pause }
};

static void vessl_e( const char * action, ... )
{
    va_list argp;
    va_start(argp, action);
    
    auto iter = vessl_actions.find(action);
    if ( iter != vessl_actions.end() )
    {
        iter->second(argp);
    }
    else
    {
        printf("[vessl] i don't know how to %s\n", action);
    }
    
    va_end(argp);
}

static bool is_prop(const char* pname, const char* name)
{
    return strcmp(pname, name)==0;
}

static vessl_value vessl_g( const char * object, const char * property )
{
    vessl_value result = { 0, VT_NONE };
    if ( is_prop(object, "output") )
    {
        if ( is_prop(property, "sr") )
        {
            if ( vessl_out.isStreamOpen() )
            {
                result.u = vessl_out.getStreamSampleRate();
                result.t = VT_UINT;
            }
        }
        else if ( is_prop(property, "bs") )
        {
            result.l = vessl_out_bufferFrames;
            result.t = VT_LONG;
        }
        else if ( is_prop(property, "t") )
        {
            if ( vessl_out.isStreamOpen() )
            {
                result.d = vessl_out.getStreamTime();
                result.t = VT_DOUBLE;
            }
        }
    }
    else
    {
        auto iter = vessl_files.find(object);
        if ( iter != vessl_files.end() )
        {
            if ( is_prop(property, "p") )
            {
                result.l = (long)sf_seek(iter->second->file, 0, SEEK_CUR);
                result.t = VT_LONG;
            }
            else if (is_prop(property, "t") )
            {
                sf_count_t seek = sf_seek(iter->second->file, 0, SEEK_CUR);
                result.d = (double)seek/(iter->second->info->samplerate);
                result.t = VT_DOUBLE;
            }
            else if ( is_prop(property, "sr") )
            {
                result.u = (unsigned int)iter->second->info->samplerate;
                result.t = VT_UINT;
            }
            else if ( is_prop(property, "ch") )
            {
                result.i = iter->second->info->channels;
                result.t = VT_INT;
            }
            else if ( is_prop(property, "f") )
            {
                result.l = (long)iter->second->info->frames;
                result.t = VT_LONG;
            }
        }
    }
    
    return result;
}

static int vessl_s( const char* object, const char* property, vessl_value value)
{
    int result = 1;
    
    if ( is_prop(object, "output") )
    {
        if ( is_prop(property, "sr") )
        {
            if ( vessl_out.isStreamOpen() == false )
            {
                vessl_out_sampleRate = value.u;
                result = 0;
            }
            else
            {
                printf("[vessl] can't set output sr because output is already open\n");
            }
        }
        else
        {
            printf("[vessl] output does not have the %s property\n", property);
        }
    }
    else if ( vessl_file* vfile = get_file(object) )
    {
        if ( is_prop(property, "p") )
        {
            if ( sf_seek(vfile->file, (sf_count_t)value.l, SEEK_SET) == -1 )
            {
                printf("[vessl] couldn't set file position because %s\n", sf_strerror(vfile->file));
            }
            else
            {
                result = 0;
            }
        }
        else if ( is_prop(property, "t") )
        {
            sf_count_t pos = (sf_count_t)(value.d * vfile->info->samplerate);
            if ( sf_seek(vfile->file, pos, SEEK_SET) == -1 )
            {
                printf("[vessl] couldn't set file time because %s\n", sf_strerror(vfile->file) );
            }
            else
            {
                result = 0;
            }
        }
        else
        {
            printf("[vessl] files do not have a %s property\n", property);
        }
    }
    else if ( vessl_expression* vexp = get_expr(object) )
    {
        vexp->eval.SetVar(property[0], value.u);
        result = 0;
    }
    
    return result;
}

std::map<std::string, std::string> vessl_defines;

void vessl_d(const char* alias, const char* value)
{
    vessl_defines[alias] = value;
}

const char* resolve_define(const char* lookup)
{
    auto iter = vessl_defines.find(lookup);
    if ( iter != vessl_defines.end() )
    {
        return iter->second.c_str();
    }
    
    return lookup;
}

vessl_struct const vessl = { vessl_e, vessl_g, vessl_s, vessl_d };
