# LCM_RNG_Crack
Break Linear Congruent Modular RNG, obtaining a, k and M from a list of output.

Simple test program at this time.  Simply type make  and then ./linmod_rev
The program will run, cracking one LCM. There is a 2nd which ATM, the code
will not crack. That 2nd one is VC's rand() function. That function is simply
(full_rand() >> 16) & 0x7fff.

This program will be a WIP, working to crack the MySpace 'auto-generated' user
accounts in the hack dump of the MySpace DB.  Those hashes have the form of
sha1($u.$p) where $p is:  FBOBH_?1?1?1?1?1?1?1?1?1?1 with ?1 being [a-yA-Y0-8]
and 10 symbols.  Each 'group'  [a-y] [A-Y] [0-8] and 10 symbols have a 25%
chance of happening for each charcater.  Then once the group is determined,
any character FROM that group is uniformly possible.  So to me, this APPEARS
to be some simplistic RNG function, and likely something like this:

char next() {
   int i = rand();
   switch(i%4) {
      case 0:  return 'a'+i%25;
      case 1:  return 'A'+i%25;
      case 2:  return '0'+i%9;
      case 3:  return symbols[i%10];
}

Now, I have no knowledge of the real sorce code which generated these 10 character
long 'random' strings. BUT it appears to be something like the function listed
above.  But this 'program', while a nice 'EXAMPLE' of how to crack a Linear
Congruent Modular pRNG (and can stand alone as that example), really is just
a stepping stone to this final goal.
