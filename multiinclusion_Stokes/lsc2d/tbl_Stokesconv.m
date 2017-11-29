% regenerate final table 5.1 for lsc2d. Barnett 10/5/14

'ext dir:', for N=100:50:350, testStokesSDevalclose(N,'e','d'); end
'int dir:', for N=100:50:350, testStokesSDevalclose(N,'i','d'); end
'ext neu:', for N=100:50:350, testStokesSDevalclose(N,'e','s'); end
'int neu:', for N=100:50:350, testStokesSDevalclose(N,'i','s'); end
