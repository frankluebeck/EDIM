/* made by 'gac', can be thrown away */
#include "src/compiled.h" 
#ifndef AVOID_PRECOMPILED
extern StructInitInfo * Init__methsel1 ( void );
extern StructInitInfo * Init__type1 ( void );
extern StructInitInfo * Init__filter1 ( void );
extern StructInitInfo * Init__oper1( void );
extern StructInitInfo * Init__random( void );
#endif
extern StructInitInfo * Init__ediv ( void );
InitInfoFunc CompInitFuncs [] = {
#ifndef AVOID_PRECOMPILED
    Init__methsel1,
    Init__type1,
    Init__oper1,
    Init__filter1,
    Init__random,
#endif
    Init__ediv,
    0
};
