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
  ; **  Calculation of the deresolution function (gaussian)     **
  ; **                                                          **
  ; **************************************************************

  function gauss_nn,ampl,ecart,moy,x

    ; |------------------------------------------------------------|
    ; |  Parameters :                                              |
    ; |  - in   moy : mean of the gaussian                         |
    ; |         ecart : standard deviation of the gaussian         |
    ; |         x : current wavelength or wavenumber               |
    ; |                                                            |
    ; |  - out  gauss : value of the factor of deresolution        |
    ; |                                                            |
    ; |------------------------------------------------------------|

    gauss_val=ampl
    gauss_val=gauss_val*exp(-ecart*(x-moy)^2)
    return, gauss_val

  end 



pro fct_deresol_block2, p, D_OPLANE, err_flag,err_msg, $
                        kresol=kresol,kunit1=kunit1,kunit2=kunit2, $
                        n_s1=n_s1, n_s2=n_s2

;This function convolves every spectrum in 'obj_spect_hr' with the CRISM response corresponding to the 'n_s' column.

  common share_dere, tab_lun, pmax, qmax, nbwp, nb_block, recpmax, wave_cm, waveref, pixa, PSF_inst, tab_pt, tab_ind 

; initialisation des tableaux
  chunk=fltarr(nbwp,qmax)
  nb_wave=long(n_elements(wave_cm))
  PSF_size=size(PSF_inst)
  nb_wref=PSF_size[2]
  nbc_block=-1
  ind_valid_q=where(pixa(*,p) ne 0.)

;INPUTS:
                                ;- 'obj_spect_hr' contains
                                ;'nbSpectra*nb_class' oversampled
                                ;spectra with 'nb_wave' each
                                ;- 'wave' array of wavenumbers
                                ;- 'n_s1:n_s2' columns of interest
                                ;- 'PSF_inst' a table containing the
                                ;    PSF parameters of the instrument
                                ; kresol type of PSF
;OUTPUTS:
                                ;- 'obj_spect_br' contains 'nbSpectra' sampled to the CRISM sampling ('nb_wref' samples) of the 'n_s' column
;CODE:

; parameters
  prec=1.e-8

; wave (cm-1 -> microns)
  wave_mi=1./wave_cm*1.e4

;Starting convolution...

  for n_wref=nb_wref-1, 0, -1  do begin

; for the first line the PSF for each wavelength are extracted
     if p eq 0 then begin

        current_PSF=replicate(0.,nb_wave)

                                ;PSF for the current band/column
                                ;couple is built...
        case kresol of
; simple gaussian PSF
           1: begin
              ind_med=fix((1./PSF_inst((n_s1+n_s2)/2,n_wref,0)*1.e4-wave_cm[0])*nb_wave/(wave_cm[nb_wave-1]-wave_cm[0]))-1
              i=1
              current_PSF(ind_med-i)=total(gauss(PSF_inst(n_s1:n_s2,n_wref,0),PSF_inst(n_s1:n_s2,n_wref,1),wave_mi(ind_med-i)))/(n_s2-n_s1+1)
              current_PSF(ind_med)=total(gauss(PSF_inst(n_s1:n_s2,n_wref,0),PSF_inst(n_s1:n_s2,n_wref,1),wave_mi(ind_med)))/(n_s2-n_s1+1)
              current_PSF(ind_med+i)=total(gauss(PSF_inst(n_s1:n_s2,n_wref,0),PSF_inst(n_s1:n_s2,n_wref,1),wave_mi(ind_med+i)))/(n_s2-n_s1+1)
              while (abs(current_PSF(ind_med-i)*current_PSF(ind_med+i)) gt prec) do begin
                 i=i+1
                 current_PSF(ind_med-i)=total(gauss(PSF_inst(n_s1:n_s2,n_wref,0),PSF_inst(n_s1:n_s2,n_wref,1),wave_mi(ind_med-i)))/(n_s2-n_s1+1)
                 current_PSF(ind_med+i)=total(gauss(PSF_inst(n_s1:n_s2,n_wref,0),PSF_inst(n_s1:n_s2,n_wref,1),wave_mi(ind_med+i)))/(n_s2-n_s1+1)
              endwhile

           end
; simple triangle PSF
           2: begin
              for s=0, nb_wave-1 do $
                 current_PSF(s)=total(triangle(PSF_inst(n_s1:n_s2,n_wref,0),PSF_inst(n_s1:n_s2,n_wref,1),wave_mi(s)))/(n_s2-n_s1+1)
           end
; simple trapezoidal PSF
           3: begin
              for s=0, nb_wave-1 do $
                 current_PSF(s)=total(trapeze(PSF_inst(n_s1:n_s2,n_wref,0),PSF_inst(n_s1:n_s2,n_wref,1),PSF_inst(n_s1:n_s2,n_wref,2),wave_mi(s)))/(n_s2-n_s1+1) 
           end
; CRISM PSF
           4: begin
              ind_med=fix((1./PSF_inst((n_s1+n_s2)/2,n_wref,3)*1.e7-wave_cm[0])*nb_wave/(wave_cm[nb_wave-1]-wave_cm[0]))-1
              i=1
              current_PSF(ind_med-i)=total(gauss_nn(PSF_inst(n_s1:n_s2,n_wref,1),PSF_inst(n_s1:n_s2,n_wref,2),PSF_inst(n_s1:n_s2,n_wref,3),wave_mi(ind_med-i)*1.e3) $
                                           +gauss_nn(PSF_inst(n_s1:n_s2,n_wref,4),PSF_inst(n_s1:n_s2,n_wref,5),PSF_inst(n_s1:n_s2,n_wref,6),wave_mi(ind_med-i)*1.e3) $
                                           +gauss_nn(PSF_inst(n_s1:n_s2,n_wref,7),PSF_inst(n_s1:n_s2,n_wref,8),PSF_inst(n_s1:n_s2,n_wref,9),wave_mi(ind_med-i)*1.e3))/(n_s2-n_s1+1)
              current_PSF(ind_med)=total(gauss_nn(PSF_inst(n_s1:n_s2,n_wref,1),PSF_inst(n_s1:n_s2,n_wref,2),PSF_inst(n_s1:n_s2,n_wref,3),wave_mi(ind_med)*1.e3) $
                                         +gauss_nn(PSF_inst(n_s1:n_s2,n_wref,4),PSF_inst(n_s1:n_s2,n_wref,5),PSF_inst(n_s1:n_s2,n_wref,6),wave_mi(ind_med)*1.e3) $
                                         +gauss_nn(PSF_inst(n_s1:n_s2,n_wref,7),PSF_inst(n_s1:n_s2,n_wref,8),PSF_inst(n_s1:n_s2,n_wref,9),wave_mi(ind_med)*1.e3))/(n_s2-n_s1+1)
              current_PSF(ind_med+i)=total(gauss_nn(PSF_inst(n_s1:n_s2,n_wref,1),PSF_inst(n_s1:n_s2,n_wref,2),PSF_inst(n_s1:n_s2,n_wref,3),wave_mi(ind_med+i)*1.e3) $
                                           +gauss_nn(PSF_inst(n_s1:n_s2,n_wref,4),PSF_inst(n_s1:n_s2,n_wref,5),PSF_inst(n_s1:n_s2,n_wref,6),wave_mi(ind_med+i)*1.e3) $
                                           +gauss_nn(PSF_inst(n_s1:n_s2,n_wref,7),PSF_inst(n_s1:n_s2,n_wref,8),PSF_inst(n_s1:n_s2,n_wref,9),wave_mi(ind_med+i)*1.e3))/(n_s2-n_s1+1)
              while (abs(current_PSF(ind_med-i)*current_PSF(ind_med+i)) gt prec) do begin 
                 i=i+1
                 current_PSF(ind_med-i)=total(gauss_nn(PSF_inst(n_s1:n_s2,n_wref,1),PSF_inst(n_s1:n_s2,n_wref,2),PSF_inst(n_s1:n_s2,n_wref,3),wave_mi(ind_med-i)*1.e3) $
                                              +gauss_nn(PSF_inst(n_s1:n_s2,n_wref,4),PSF_inst(n_s1:n_s2,n_wref,5),PSF_inst(n_s1:n_s2,n_wref,6),wave_mi(ind_med-i)*1.e3) $
                                              +gauss_nn(PSF_inst(n_s1:n_s2,n_wref,7),PSF_inst(n_s1:n_s2,n_wref,8),PSF_inst(n_s1:n_s2,n_wref,9),wave_mi(ind_med-i)*1.e3))/(n_s2-n_s1+1)
                 current_PSF(ind_med+i)=total(gauss_nn(PSF_inst(n_s1:n_s2,n_wref,1),PSF_inst(n_s1:n_s2,n_wref,2),PSF_inst(n_s1:n_s2,n_wref,3),wave_mi(ind_med+i)*1.e3) $
                                              +gauss_nn(PSF_inst(n_s1:n_s2,n_wref,4),PSF_inst(n_s1:n_s2,n_wref,5),PSF_inst(n_s1:n_s2,n_wref,6),wave_mi(ind_med+i)*1.e3) $
                                              +gauss_nn(PSF_inst(n_s1:n_s2,n_wref,7),PSF_inst(n_s1:n_s2,n_wref,8),PSF_inst(n_s1:n_s2,n_wref,9),wave_mi(ind_med+i)*1.e3))/(n_s2-n_s1+1)
              endwhile
           end
        endcase
; extraction of the spectral range where PSF admits non zero values 
        ext_PSF=current_PSF[ind_med-i:ind_med+i]
; keep that into memory for future uses across procedure calls
        tab_ind[n_wref,0]=ind_med-i
        tab_ind[n_wref,1]=ind_med+i
        tab_pt[n_wref]=ptr_new(ext_PSF)
        kmin=ind_med-i
        kmax=ind_med+i
     endif else begin
        kmin=tab_ind[n_wref,0]
        kmax=tab_ind[n_wref,1]
        ext_PSF=(*tab_pt[n_wref])
     endelse

     ext_spect=fltarr(kmax-kmin+1,qmax)

                                ;the convolution with all spectra in
                                ;'obj_spect_hr' is performed by the
                                ;dot product

     for k=kmin,kmax do begin

        nbc_blockp=nbc_block
                                ; calcul du no de bloc correspondant à la longueur d'onde courante
        nbc_block=fix((float(k+1)-0.1)/float(nbwp))
                                ; chargement en mémoire du bloc courant s'il ne l'est pas déjà
        if(nbc_block ne nbc_blockp) then begin
                                ; calcul du no de fichier scratch pour le bloc courant
           sf_nb=(p*nb_block+nbc_block)/recpmax
           rec_nb=((p*nb_block+nbc_block) MOD recpmax)
           data=assoc(tab_lun[sf_nb],chunk)
           pix_mem=data(rec_nb)
;		print,' Nb block: ',nbc_block
        endif

                                ; calcul de l'indice correspondant à la longueur d'onde courante pour le bloc courant
        nbcwp=(k MOD nbwp)
; extraction of the spectral range where PSF admits non zero values
        ext_spect[k-kmin,ind_valid_q]=pix_mem(nbcwp,ind_valid_q)/pixa(ind_valid_q,p)

     endfor  

; eliminating NaN values (if any) for the calculation
     for ind=0,n_elements(ind_valid_q)-1 do begin
        q=ind_valid_q[ind]
        ind_valid_w=where(finite(ext_spect[*,q]))
        response_inst=total(replicate(1,n_elements(ind_valid_w))*ext_PSF[ind_valid_w])
        D_OPLANE(n_wref,q)=total(ext_spect[ind_valid_w,q]*ext_PSF[ind_valid_w])/$
                           response_inst
     endfor

  endfor  

  err_flag=0

  return

end
