
      // // -- upwind scheme-- 
      // convec_scheme = 1;
      // aE[i][j] = De + fmax(-Fe, 0.0);
      // aW[i][j] = Dw + fmax(Fw, 0.0);
      // aN[i][j] = Dn + fmax(-Fn, 0.0);
      // aS[i][j] = Ds + fmax(Fs, 0.0);
     
      // -- C-D scheme--
      aE[i][j] = De - 0.5 * Fe;
      aW[i][j] = Dw + 0.5 * Fw;
      aN[i][j] = Dn - 0.5 * Fn;
      aS[i][j] = Ds + 0.5 * Fs;

      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j] - Sp[i][j] + Fe + Fn - Fw - Fs;
