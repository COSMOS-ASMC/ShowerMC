BEGIN {nl = 0; nv=0; ntv=0; nfc=0; first=0;
      filec=0; prev=" "; prevline=" "; contline=" ";
      system("echo ' ' > tempvectvtx");
      system("echo ' ' > tempvectnv");
      system("echo ' ' > tempvectcolor");
 }
#  skip until 2nd NF==0 
nfc < 2 && NF == 0 {nfc++;next};
nfc < 2 {next};

NF==3 {ntv++; nv++; prev="vtx"; prevline=$0; print $0 >> "tempvectvtx"}
ntv >= maxsize &&  nv>1  { if(prev != "sep")
			     { nl++;
			       system("echo ' ' >> tempvectvtx");
			       print nv >> "tempvectnv"; 
			       if(first == 1)   print "0" >> "tempvectcolor";
			       else  print "1" >>  "tempvectcolor";
			     }
       filec++;
       filex = file filec".vect";					 
       print "VECT" > filex;
       print  nl, ntv, "1" >> filex;
       system("cat tempvectnv >>" filex);
       system("cat tempvectcolor >>" filex);


       if( filec  > 1 && contline != " " ) {    
	 system("echo " contline ">>" filex)
       }

       contline = prevline;

   			   
       system("cat tempvectvtx >>" filex);
       system("echo ' ' >>" filex) ;


       system("echo "  color  ">>" filex); 

       ntv=1; nv=1; nl=0; first=0;

      system("echo ' ' > tempvectvtx");
      system("echo ' ' > tempvectnv");
      system("echo ' ' > tempvectcolor");
       }

NF==0 && prev=="vtx" {prev="sep"; if(ntv == 1)
			{ contline = " "; nv--; ntv--}
                        else { nl++;
			       system("echo ' ' >> tempvectvtx");
			       print nv >> "tempvectnv"; 
			       if(first == 1)   print "0" >> "tempvectcolor";
			       else  print "1" >>  "tempvectcolor";
			       first=1; nv=0; 
			     }
		    }

END {

  if( ntv > 0 ) {
    filec++;
    filex = file filec".vect";
    print "VECT" > filex;
    print  nl, ntv, "1" >> filex;

    system("cat tempvectnv >>" filex);
    system("echo ' ' >>" filex) ;
    system("cat tempvectcolor >>" filex);
    system("echo ' ' >>" filex) ;

    
    if( filec  > 1 && contline != " " ) {    
      system("echo " contline ">>" filex)
    }

    system("cat tempvectvtx >>" filex);
    system("echo ' ' >>" filex) ;
    system("echo "  color  ">>" filex); 

    system("rm -f tempvectvtx");				 
    system("rm -f tempvectnv");				 
    system("rm -f tempvectcolor");
  }
}
 

