#ifndef CEC_NORET_H
#define CEC_NORET_H
#if __STDC_VERSION__ >= 201112L
    #include <stdnoreturn.h>
#elif __GNUC__ >= 3
    #define noreturn
#endif
#endif //CEC_NORET_H
