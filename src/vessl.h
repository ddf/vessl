//
//  vessl.h
//  vessl
//
//  Created by Damien Di Fede on 10/19/15.
//  Copyright © 2015 compartmental. All rights reserved.
//

#ifndef vessl_h
#define vessl_h

typedef int (* vessl_render_callback)(float* buffer, unsigned int nBufferFrames, void* userData);
typedef int (* vessl_render_unit)(float* sampleFrame, void* userData);

typedef enum
{
    VT_NONE,
    VT_INT,
    VT_UINT,
    VT_LONG,
    VT_FLOAT,
    VT_DOUBLE
} vessl_type;

typedef struct
{
    union
    {
        int i;
        unsigned int u;
        long l;
        float f;
        double d;
    };
    
    vessl_type t;
} vessl_value;

typedef struct
{
    // execute an action
    void (* const e)(const char * action, ...);
    
    // get the value of a property, t will be VT_NONE if this fails
    vessl_value (* const g)(const char * object, const char * property);

    // set the value of a property, returns 0 if it worked, 1 if not
    int (* const s)(const char * object, const char * property, vessl_value value);
    
    // define an alias, works like #define
    void (* const d)(const char* alias, const char* value);
} vessl_struct;

extern vessl_struct const vessl;


#endif /* vessl_h */
