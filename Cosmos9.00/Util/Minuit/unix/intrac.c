/* int intrac_() */
/*  absoft @ mac */
#if defined (MACOSX)
 int INTRAC() 
#else
/* int INTRAC_()   */
 int intrac_()  
#endif

{
    return ((int) isatty(0));
}
