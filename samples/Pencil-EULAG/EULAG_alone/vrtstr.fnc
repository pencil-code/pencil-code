!     fs(zt) = istr*(exp((zb - zt)*log(1. + SD)/zb) - 1.)/SD
      fs(za)=  istr*(1.-exp(-za/SD))             +(1-istr)*za
      fsi(zt)=-istr*SD*alog(1.-istr*zt/zb*fs(zb))+(1-istr)*zt
      fsk(zk)=zb/fs(zb)*fs(zk)
