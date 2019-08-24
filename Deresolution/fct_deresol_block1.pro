  ; **************************************************************
  ; **                                                          **
  ; **  Calculation of the deresolution function (gaussian)     **
  ; **                                                          **
  ; **************************************************************

  function gauss,moy,ecart,x

    ; |------------------------------------------------------------|
    ; |  Parameters :                                              |
    ; |  - in   moy : mean of the gaussian                         |
    ; |         ecart : standard deviation of the gaussian         |
    ; |         x : current wavelength or wavenumber               |
    ; |                                                            |
    ; |  - out  gauss : value of the factor of deresolution        |
    ; |                                                            |
    ; |------------------------------------------------------------|

    gauss_val=1./(sqrt(2*!pi)*ecart)
    gauss_val=gauss_val*exp(-(x-moy)^2/(2*ecart^2))
    return, gauss_val

  end 

  ; **************************************************************
  ; **                                                          **
  ; **  Calculation of the deresolution function (triangle)     **
  ; **                                                          **
  ; **************************************************************

  function triangle,moy,ecart,x

    ; |------------------------------------------------------------|
    ; |  Parameters :                                              |
    ; |  - in   moy : mean of the triangle                         |
    ; |         ecart : standard deviation of the triangle         |
    ; |         x : current wavelength or wavenumber               |
    ; |                                                            |
    ; |  - out  triangle : value of the factor of deresolution     |
    ; |                                                            |
    ; |------------------------------------------------------------|


    if(abs(x-moy) GT (2.*ecart)) then triangle_val=0.  else $
   if((x-moy) LE 0.) then triangle_val=(x-moy)+2.*ecart else $ 
      triangle_val=(moy-x)+2.*ecart
     triangle_val=triangle_val/(4.*ecart^2)

    return, triangle_val
  end 


  ; **************************************************************
  ; **                                                          **
  ; **  Calculation of the deresolution function (trapeze)      **
  ; **                                                          **
  ; **************************************************************

  function trapeze,moy,ecart,top,x

    ; |------------------------------------------------------------|
    ; |  Parameters :                                              |
    ; |  - in   moy : mean of the trapeze                          |
    ; |         ecart : standard deviation of the trapeze          |
    ; |         top : half width at the top of the trapeze         |
    ; |         x : current wavelength or wavenumber               |
    ; |                                                            |
    ; |  - out  trapeze : value of the factor of deresolution      |
    ; |                                                            |
    ; |------------------------------------------------------------|

    if(abs(x-moy) GT (2.*ecart-top)) then trapeze_val=0.$
    else begin
       if(abs(x-moy) LE top) then trapeze_val=1/(2*ecart)$
       else begin
          if((x-moy) LE 0.) then trapeze_val=((x-moy)+2.*ecart-top)/(4.*ecart*(ecart-top))$
          else trapeze_val=((moy-x)+2.*ecart-top)/(4.*ecart*(ecart-top))
      endelse 
    endelse  

    return, trapeze_val 
  end 

  ; **************************************************************

pro fct_deresol_block1,p,D_OPLANE,err_flag,err_msg, $
                      resol=resol,topresol=topresol,kresol=kresol,kunit1=kunit1,kunit2=kunit2

common share_dere, tab_lun, pmax, qmax, nbwp, nb_block, recpmax, wavep, waverefp, pixa

; initialisation des tableaux
  chunk=fltarr(nbwp,qmax)
  nref=long(n_elements(waverefp))
  nfich=long(n_elements(wavep))
  waveref=fltarr(nref+1)
  wave=fltarr(nfich+1)

; Emulation of Fortran's norms for indexes
  waveref(1:nref)=waverefp(0:nref-1)
  wave(1:nfich)=wavep(0:nfich-1)

                                ;+
                                ;    $Id: fct_deresol_block.f90,v 1.1 2004/10/01 09:10:53 madec Exp doute $
                                ;
                                ; NAME: fct_deresol_block
                                ;
                                ;
                                ;
                                ; PURPOSE: 
                                ;	convolue un spectre brut par une fonction instrumentale
                                ; 	caract�ris�e par sa forme (gaussienne, triangle ou trap�ze) et sa
                                ; 	r�solution spectrale. Le calcul est effectu� suivant un
                                ; 	�chantillonnage � fournir (liste de longueurs d'onde).
                                ;
                                ;
                                ;
                                ; CATEGORY: 
                                ;	traitement du signal
                                ;
                                ;
                                ;
                                ; CALLING SEQUENCE: 
                                ;	call fct_deresol(nref,waveref,resol,topresol,$
                                ;                        $kresol,nfich,wave,refh,drefh,kunit1,kunit2) 
                                ;
                                ;
                                ;
                                ; INPUTS:
                                ;
                                ;   nref dimensionnalit� de l'�chantillonnage du spectre convolu�
                                ;   waveref �chantillonnage du spectre convolu� 
                                ;   nfich dimensionnalit� du spectre brut
                                ;   wave �chantillonnage du spectre brut
                                ;   refh spectre brut
                                ;
                                ;
                                ;
                                ; OPTIONAL INPUTS:
                                ;
                                ;   resol r�solution spectrale du spectre convolu� (si param�tre
                                ;     absent, la r�solution NIMS est choisie par d�faut)
                                ;   topresol r�solution spectrale du spectre convolu� (cas du
                                ;     trap�ze) (si param�tre absent, la r�solution est prise � 0.)
                                ;   kunit1, kunit2 unit�s spectrales respectives du spectre brut et
                                ;     convolu�
                                ;     kunit: cm-1 -> 1 , micron -> 2 , nm -> 3 
                                ;   kresol type of the spectral response function : 
                                ;     kresol: Gaussienne -> 1,  Triangle -> 2,  Trapeze -> 3
                                ;
                                ;
                                ; KEYWORD PARAMETERS:
                                ;
                                ;
                                ;
                                ; OUTPUTS:
                                ;
                                ;   D_OPLANE current spatial-spectral plan
                                ;   err_flag flag erreur (=0-> pas d'erreur, =1-> erreur)
                                ;   err_msg message erreur (vide si err_flag=0)
                                ;
                                ;
                                ; OPTIONAL OUTPUTS:
                                ;
                                ;
                                ;
                                ; COMMON BLOCKS:
                                ;
                                ;
                                ;
                                ; SIDE EFFECTS:
                                ;
                                ;
                                ;
                                ; RESTRICTIONS:
                                ;
                                ;
                                ;
                                ; PROCEDURE:
                                ;
                                ;
                                ;
                                ; EXAMPLE:
                                ;
                                ;
                                ;
                                ; MODIFICATION HISTORY:
                                ; 
                                ;     $Log: fct_deresol_block.f90,v $
                                ;     Revision 1.1  2004/10/01 09:10:53  madec
                                ;     Initial revision
                                ;
                                ;
                                ;              Bernard Schmitt - LPG Grenoble - Version initiale
                                ;   16/08/2001 Sylvain Dout�   - LPG Grenoble - Version f90
                                ;
                                ;-

                                ; initialisation des variables
  inb=0 & nb=0 & ind=0 & j=0 & i=0 
  q=0 & nbc_block=0 & nbc_blockp=0 & qpmax=0 & sf_nb=0 & nbcwp=0 & rec_nb=0
  x=0. & sigma=0. & topsigma=0. & pas=0. & wavetrinf=0. & wavetrsup=0. & waveinf=0. & wavesup=0.
  xinf=0. & xsup=0. & wavelm=0. & wavemm=0. & wavel=0. & wavemp=0. & wavelp=0. & drefi=0. & waveb=0.
  wavebm=0. & wavebp=0. & refj=0.
  pasnu=0. & wavenuref=0. & sigmanu=0. 
                                ; parameter(x=5.)
  nbf=50000L
  err_msg1=' ' & err_msg2=' ' & err_msg3=' ' & err_msg4=' '
  err_msg5=' ' & err_msg6=' ' & err_msg7=' '

                                ; *****************************************************
                                ; Initialisation du numero de bloc courant

  nbc_block=0

                                ; *****************************************************

                                ; Traitement des param�tres facultatifs

  if(keyword_set(kresol) eq 0) then kresol=2
  if(keyword_set(resol) eq 0) then resol=0.024
  if(keyword_set(topresol) eq 0) then topresol=0.
  if(keyword_set(kunit1) eq 0) then kunit1=1
  if(keyword_set(kunit2) eq 0) then kunit2=2

  if (kunit1 EQ kunit2) then begin

                                ; nb de point echantillonnant la fonction de deresolution 
                                ; (x: facteur de coupure)
                                ; *******************************************************

     sigma=resol/2.
     topsigma=topresol/2.
     pas=wave(2)-wave(1)
     case kresol of
        1: x=5.*(1+pas/resol)
        2: x=2.
        else: x=2.-topresol/resol
     endcase
     inb=fix(x*sigma/pas)
     nb=2*inb+1
                                ; print,'nb of calculation points : ',nb
                                ; print,'Cut factor: x = ',x
     ind=1

     for i=1,nref do begin
                                ; ****************************************************
                                ; Limites du triangle d'integration
                                ; ****************************************************
        wavetrinf=waveref(i)-2.*sigma
        wavetrsup=waveref(i)+2.*sigma

                                ; recherche de la wave juste superieure a waveref
                                ; *****************************************************
        for j=ind,nfich do begin
           if (wave(j) LE waveref(i)) then begin
              continue
           endif else begin
              ind=j
              waveinf=wave(ind-inb)
              wavesup=wave(ind+inb)
              if ((waveinf LT wave(1)) or (wavesup GT wave(nfich))) then begin
                 xinf=(wave(1)-waveinf)/sigma
                 xsup=(wavesup-wave(nfich))/sigma
                 err_msg1=''
                 err_msg1='current value of the Gauss factor x: '+string(x)
                 err_msg2='current inf. limit : '+string(wave(1))
                 err_msg3='needed inf. limit : '+string(waveinf)
                 err_msg4='needed value for x(inf) : '+string(xinf)
                 err_msg5='current sup. limit : '+string(wave(nfich))
                 err_msg6='needed sup. limit : '+string(wavesup)
                 err_msg7='needed value for x(sup) : '+string(xsup)
                 err_msg='Spectral range of the file to be convolved too small '+$
                         strcompress(err_msg1)+strcompress(err_msg2)+strcompress(err_msg3)+$
                         strcompress(err_msg4)+strcompress(err_msg5)+strcompress(err_msg6)+$
                         strcompress(err_msg7) 
                 err_flag=1
                 return
              endif 
              break
           endelse  
        endfor  

                                ; convolution fonction deresolution entre waveref +/- x*sigma
                                ; ***********************************************************

        D_OPLANE(i-1,*)=0.

        for j=ind-inb-1,ind+inb+1 do begin

           nbc_blockp=nbc_block
                                ; calcul du no de bloc correspondant � la longueur d'onde courante
           nbc_block=fix((float(j)-0.1)/float(nbwp))+1
                                ; chargement en m�moire du bloc courant s'il ne l'est pas d�j�
           if(nbc_block ne nbc_blockp) then begin
                                ; calcul du no de fichier scratch pour le bloc courant
              sf_nb=(p*nb_block+nbc_block-1)/recpmax
              rec_nb=((p*nb_block+nbc_block-1) MOD recpmax)
              data=assoc(tab_lun[sf_nb],chunk)
              pix_mem=data(rec_nb)
           endif
                                ; calcul de l'indice correspondant � la longueur d'onde courante pour le bloc courant
           nbcwp=((j-1) MOD nbwp)

           wavelm=wave(j-1)
           wavemm=(wave(j-1)+wave(j))/2
           wavel=wave(j)
           wavemp=(wave(j)+wave(j+1))/2
           wavelp=wave(j+1)

           ind_valid=where(pixa(*,p) ne 0.)
           refj=pix_mem(nbcwp,ind_valid)/pixa(ind_valid,p)


           case kresol of
              1: begin
                 drefi=refj*gauss(waveref(i),sigma,wavel)
                 D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*(wavelp-wavelm)/2
              end
              2: begin
                 case 1B of
                                ; traitement bord inferieur triangle
                                ; **********************************
                    ((wavel GE wavetrinf) and (wavelm LT wavetrinf)) : begin
                       waveb=((wavel+wavelp)/2.+wavetrinf)/2.
                       drefi=refj*triangle(waveref(i),sigma,waveb)
                       D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs((wavel+wavelp)/2-wavetrinf)
                    end
                                ; traitement bord superieur triangle
                                ; **********************************
                    ((wavel LE wavetrsup) and (wavelp GT wavetrsup)) : begin
                       waveb=((wavel+wavelm)/2.+wavetrsup)/2.
                       drefi=refj*triangle(waveref(i),sigma,waveb)
                       D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs(wavetrsup-(wavel+wavelm)/2)
                    end
                                ; traitement sommet triangle
                                ; **********************************
                    ((wavel LT waveref(i)) and (wavemp GT waveref(i))) : begin
                       wavebm=((wavel+wavelm)/2.+waveref(i))/2.
                       wavebp=((wavel+wavelp)/2.+waveref(i))/2.
                       drefi=refj*triangle(waveref(i),sigma,wavebm)
                       D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs(waveref(i)-(wavel+wavelm)/2.)
                       drefi=refj*triangle(waveref(i),sigma,wavebp)
                       D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs((wavel+wavelp)/2.-waveref(i))
                    end
                    else: begin
                                ; traitement genral triangle
                                ; **********************************
                       drefi=refj*triangle(waveref(i),sigma,wavel)
                       D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs(wavelp-wavelm)/2
                    end
                 endcase
              end  
              else: begin
                 drefi=refj*trapeze(waveref(i),sigma,topsigma,wavel)
                 D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs(wavelp-wavelm)/2
              end 
           endcase

        endfor

     endfor

  endif else begin

                                ; ****************************************************
                                ; cas du changement d'unite
                                ; ****************************************************

                                ; nb de point echantillonnant la fonction de deresolution 
                                ; (x: facteur de coupure)
                                ; *******************************************************

     sigma=resol/2.
     pasnu=wave(2)-wave(1)
     case kresol of
        1: x=5.
        2: x=2.
        else: x=2.-topresol/resol
     endcase

     for i=1,nref do begin
        wavenuref=1.e4/waveref(i)
        sigmanu=1.e4*sigma/waveref(i)^2
        inb=fix(x*sigmanu/pasnu)
        nb=2*inb+1
        wavetrinf=waveref(i)-2*sigma
        wavetrsup=waveref(i)+2*sigma

                                ; recherche de la wave juste superieure a waveref
                                ; ***********************************************
        for j=1,nfich do begin
           if (wave(j) LE wavenuref) then begin
              continue
           endif else begin
              ind=j
              case 1B of
                 (ind LE inb) :begin
                    waveinf=wavenuref-x*sigmanu
                    xinf=(wavenuref-wave(1))/sigmanu
                    err_msg1='current value of the Gauss factor x: '+string(x)
                    err_msg2='current inf. limit : '+string(wave(1))
                    err_msg3='needed inf. limit : '+string(waveinf)
                    err_msg4='needed value for x(inf) : ' +string(xinf)
                    err_msg='Spectral range of the file to be convolved too small %'$
                            +strcompress(err_msg1)+strcompress(err_msg2)+strcompress(err_msg3)$
                            +strcompress(err_msg4)
                    err_flag=1
                    return 
                 end
                 (ind GT (nfich-inb)) : begin
                    wavesup=wavenuref+x*sigmanu
                    xsup=(wave(nfich)-wavenuref)/sigmanu
                    err_msg1=''
                    err_msg1='current value of the Gauss factor x: '+string(x)
                    err_msg5='current sup. limit : '+string(wave(nfich))
                    err_msg6='needed sup. limit : '+string(wavesup)
                    err_msg7='needed value for x(sup) : '+string(xsup)
                    err_msg='Spectral range of the file to be convolved too small %'+$
                            strcompress(err_msg1)+strcompress(err_msg5)+strcompress(err_msg6)+$
                            strcompress(err_msg7)
                    err_flag=1
                    return
                 end
                 else : begin
                    waveinf=wave(ind-inb)
                    wavesup=wave(ind+inb)
                                ;   print,'inf. limit file : ',wave(1)
                                ;   print,'used inf. limit : ',waveinf
                                ;   print,'used sup. limit : ',wavesup
                                ;   print,'sup. limit file : ',wave(nfich)
                 end
              endcase
              break
           endelse
        endfor

                                ; convolution fct deresolution entre waveref +/- 1e4/x*sigmanu
                                ; ************************************************************

        D_OPLANE(i-1,*)=0.

        for j=ind-inb,ind+inb do begin

           nbc_blockp=nbc_block
                                ; calcul du no de bloc correspondant � la longueur d'onde courante
           nbc_block=fix((float(j)-0.1)/float(nbwp))+1
                                ; chargement en m�moire du bloc courant s'il ne l'est pas d�j�
           if(nbc_block ne nbc_blockp) then begin
                                ; calcul du no de fichier scratch pour le bloc courant
              sf_nb=(p*nb_block+nbc_block-1)/recpmax
              rec_nb=((p*nb_block+nbc_block-1) MOD recpmax)
              data=assoc(tab_lun[sf_nb],chunk)
              pix_mem=data(rec_nb)
              
           endif
                                ; calcul de l'indice correspondant � la longueur d'onde courante pour le bloc courant
           nbcwp=((j-1) MOD nbwp)

           wavelm=1.e4/wave(j-1)
           wavel=1.e4/wave(j)
           wavelp=1.e4/wave(j+1)


           ind_valid=where(pixa(*,p) ne 0.)
           refj=pix_mem(nbcwp,ind_valid)/pixa(ind_valid,p)

           case kresol of
              1: begin
                 drefi=refj*gauss(waveref(i),sigma,wavel)
                 D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs(wavelp-wavelm)/2
              end
              2: begin
                 case 1B of
                    ((wavel GE wavetrinf) and (wavelm LT wavetrinf)) : begin
                       waveb=((wavel+wavelp)/2.+wavetrinf)/2.
                       drefi=refj*triangle(waveref(i),sigma,waveb)
                       D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs((wavel+wavelp)/2-wavetrinf)
                    end
                    ((wavel LE wavetrsup) and (wavelp GT wavetrsup)): begin
                       waveb=((wavel+wavelm)/2.+wavetrsup)/2.
                       drefi=refj*triangle(waveref(i),sigma,waveb)
                       D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs(wavetrsup-(wavel+wavelm)/2)
                    end
                    else: begin 
                       drefi=refj*triangle(waveref(i),sigma,wavel)
                       D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs(wavelp-wavelm)/2
                    end

                 endcase
              end
              else: begin
                 drefi=refj*trapeze(waveref(i),sigma,topsigma,wavel)
                 D_OPLANE(i-1,ind_valid)=D_OPLANE(i-1,ind_valid)+drefi*abs(wavelp-wavelm)/2
              end

           endcase

        endfor 

     endfor 

  endelse 

  err_flag=0

  return 

end 





