PRO convol_inst_spectribase_block, SUFFIX, N_ENDM, C=C, L=L, keep=keep, INST_OMEGA=INST_OMEGA, INST_CRISM=INST_CRISM

  common share_dere, tab_lun, pmax, qmax, nbwp, nb_block, recpmax, tab_w, tab_wref, pixa, PSF_inst, tab_pt, tab_ind 

workdir='/data/fschmidt/fandrieu/fastworkdir/'
env='/home/fandrieu/Inversion/job/env.txt'
; *** initialization
  nb_block=0
  nbwp=0L
  nb_sfiles=0
  recpmax=0L
  nbw=0L

  err_msg=' '
;***Reading envfich
  OPENR, 1, env, ERROR=ios
  if (ios ne 0) then goto, mistake9000
  namedir=''
  READF, 1, FORMAT='(28X,A90)', namedir
  homedir=''
  READF, 1, FORMAT='(17x,A90)', homedir
  listareadir=''
  READF, 1, FORMAT='(21x,A90)', listareadir
  scenareadir=''
  READF, 1, FORMAT='(21X,A90)', scenareadir
  charactareadir=''
  READF, 1, FORMAT='(24X,A90)', charactareadir
  geoareadir=''
  READF, 1, FORMAT='(20X,A90)', geoareadir
  geopixdir=''
  READF, 1, FORMAT='(19X,A90)', geopixdir
  geovisildir=''
  READF, 1, FORMAT='(21X,A90)', geovisildir
  atmtaudir=''
  READF, 1, FORMAT='(21X,A90)', atmtaudir
  charaspectrdir=''
  READF, 1, FORMAT='(24X,A90)', charaspectrdir
  derefichdir=''
  READF, 1, FORMAT='(30x,A90)', derefichdir
  close, 1

;***Opening the spectral characteristics file
  case 1B of
    (keyword_set(C) eq 1): CHANNEL='C'
    (keyword_set(L) eq 1): CHANNEL='L'
    else:  CHANNEL='V'
  endcase

  charaspectr='spect_' + strcompress(SUFFIX, /remove_all) + '_em' + strcompress(N_ENDM, /remove_all) + '_' + CHANNEL + '.txt'

  charaspectr = strcompress(homedir, /remove_all) + strcompress(charaspectrdir, /remove_all) + '/' + charaspectr

  OPENR, 2, charaspectr, ERROR=ios
  if (ios ne 0) then goto, mistake9095

;*****Opening the raw binary scratch files

  scratch_file = '_'+strcompress(SUFFIX, /remove_all) + '_em' + strcompress(N_ENDM, /remove_all) + '_' + CHANNEL + '_' + 'scratch_header'
  scratch_file = workdir +  scratch_file

  if(keyword_set(keep) eq 1) then OPENR, 3, scratch_file, ERROR=ios, /F77_UNFORMATTED else OPENR, 3, scratch_file, ERROR=ios, /F77_UNFORMATTED, /delete
  if (ios ne 0) then goto, mistake9096

  readu, 3, nb_block,nbwp,nb_sfiles,recpmax,nbw
  tab_w = fltarr(nbw)           ;tab_w is the raw list of wavenumbers
  refhp=fltarr(nbw)
  readu, 3, tab_w

  tab_lun=intarr(nb_sfiles)
  for n_sfiles=0,nb_sfiles-1 do begin
     scratch_file = '_'+strcompress(SUFFIX, /remove_all) + '_em' + strcompress(N_ENDM, /remove_all) + '_' + CHANNEL + '_' + 'scratch'+$
                    strcompress(n_sfiles+1, /remove_all)
     scratch_file = workdir + scratch_file
     if(keyword_set(keep) eq 1) then openr, lun, scratch_file, /get_lun else openr, lun, scratch_file, /delete, /get_lun  ;only 28 lun available
     tab_lun[n_sfiles]=lun 
;     lun=nb_sfiles+n_sfiles   ;to be removed
;     if(keyword_set(keep) eq 1) then openr, lun, scratch_file else openr, lun, scratch_file, /delete   ;to be removed
;     tab_lun[n_sfiles]=lun    ;to be removed
  endfor

  final_output = strcompress(SUFFIX, /remove_all) + '_em' + strcompress(N_ENDM, /remove_all) + '_' + CHANNEL + '.cub'
  final_output = strcompress(homedir, /remove_all) + strcompress(derefichdir, /remove_all) + '/' + final_output

;;   raw_output = strcompress(SUFFIX, /remove_all) + '_em' + strcompress(N_ENDM, /remove_all) + '_' + CHANNEL + '.output.raw'
;;   raw_output = strcompress(atmtaudir, /remove_all) + '/' + raw_output
;;   final_output = strcompress(SUFFIX, /remove_all) + '_em' + strcompress(N_ENDM, /remove_all) + '_' + CHANNEL
;;   final_output = strcompress(atmtaudir, /remove_all) + '/' + final_output

;***Reading the parameters in the spectral characteristics file
  dim1=''
  readf, 2, FORMAT='(54x,A2)', dim1

  case dim1 of
     'mi': unit1=2
     'cm': unit1=1
     'nm': unit1=3
     else: goto, mistake9100
  endcase

  if (dim1 eq 'mi') then begin
     readf, 2, FORMAT='(48x,f7.4,1x,f7.4)', fww, lww
  endif else begin
     readf, 2, FORMAT='(48x,f7.1,1x,f7.1)', fww, lww
  endelse

  dim2='dim2'
  readf, 2, FORMAT='(52x,A2)', dim2

  case dim2 of
     'mi': begin
        unit2=2
        str_unit2='MICROMETER'
        str_length2=10
     end
     'cm': begin
        unit2=1
        str_unit2='WAVENUMBER'
        str_length2=10
     end
     'nm': begin
        unit2=3
        str_unit2='NANOMETER'
        str_length2=9
     end
     else: goto, mistake9105
  endcase

  readf, 2, FORMAT='(42x,I5)', pmax
  print, FORMAT='("Nb of image lines : ",I5)', pmax
  readf, 2, FORMAT='(44x,I5)', qmax
  print, FORMAT='("Nb of image samples : ",I5)', qmax

  pixa=dblarr(qmax,pmax)
  readu, 3, pixa

;;   * image option
  opt3='opt3'
  readf, 2, FORMAT='(15x,A4)', opt3
  if ((opt3 ne 'raw') and (opt3 ne 'conv')) then goto, mistake9111

;;   * correction option
  opt4='opt4'
  readf, 2, FORMAT='(20x,A1)', opt4
  if ((opt4 ne 'n') and (opt4 ne 's') and (opt4 ne 'd')) then goto, mistake9112

;;   * interpolation option
  opt5='opt5'
  readf, 2, FORMAT='(23x,A1)', opt5
  if ((opt5 ne 'y') and (opt5 ne 'n')) then goto, mistake9113

;;   * number of Fourier components used
  readf, 2, FORMAT='(36x,I3)', nbord
  expo=fltarr(nbord+1)

  readf, 2, FORMAT='(37x,I1)', instfct
  print, FORMAT='("Type of instrument function : ",I1)', instfct

  case 1B of
    (keyword_set(INST_CRISM) eq 1): instfct=4
    (keyword_set(INST_OMEGA) eq 1): instfct=instfct;2
    else:  instfct=4
  endcase
  

  case 1B of
     ;(instfct eq 2): begin MODIF F.ANDRIEU 09/2014
     (instfct NE 4): begin
        readf, 2, FORMAT='(22x,f9.5)', dlamb
        print, FORMAT='("  spectral resolution : ",f9.5)', dlamb
        if (instfct eq 3) then begin
           readf, 2, FORMAT='(26x,f9.5)', topdlamb
           print, FORMAT='("  top spectral resolution : ",f9.5)', topdlamb
        endif else begin
           topdlamb=0.
        endelse
        derefich=''
        readf, 2, FORMAT='(20x,A28)', derefich
        print, FORMAT='("  deresolution file : ",A28)', derefich
        if(derefich eq 'none') then begin
           readf, 2, FORMAT='(20x,f9.5)', slamb
           readf, 2, FORMAT='(33x,E11.5)', firstw
           readf, 2, FORMAT='(32x,E11.5)', lastw
           print, FORMAT='("  spectral sampling : ",f9.5)', slamb
           print, FORMAT='("  first wavelength : ",E10.5)', firstw
           print, FORMAT='("  last wavelength : ",E10.5)', lastw
           nbwref=floor(abs(lastw-firstw)/slamb)+1
           tab_wref = fltarr(nbwref)
           tab_drefl = fltarr(nbwref)
           tab_dlamb = fltarr(nbwref)
           tab_index = fltarr(nbwref)
           for nbcw=0, nbwref-1 do begin
              tab_wref(nbcw)=(lastw-firstw)/(nbwref-1)*(nbcw-1)+firstw
              tab_dlamb(nbcw)=dlamb
           endfor
           close, 2
        endif else begin
           derefich = strcompress(homedir, /remove_all) + strcompress(derefichdir, /remove_all) + '/' + strcompress(derefich, /remove_all)
           openr, 4, derefich, ERROR=ios
           if (ios ne 0) then goto, mistake9150
           ;; READF, 4, FORMAT='(f7.4)', nbwref
           READF, 4, nbwref
           tab_wref = fltarr(nbwref)
           tab_drefl = fltarr(nbwref)
           tab_dlamb = fltarr(nbwref)
           tab_index = fltarr(nbwref)
           for nbcw=0, nbwref-1 do begin
              readf, 4, wave
              tab_wref(nbcw)=wave
              tab_dlamb(nbcw)=dlamb
           endfor
           close, 2, 4
        endelse
        output_cube = fltarr(nbwref,qmax,pmax)
        D_OPLANE=fltarr(nbwref,qmax)
        for p=0,pmax-1 do begin
           print,format='("line ",I5," treated")',p
           fct_deresol_block1,p,D_OPLANE,err_flag,err_msg,resol=dlamb,topresol=topdlamb,kresol=instfct,kunit1=unit1,kunit2=unit2
           if(err_flag ne 0) then begin
              goto, mistake9170
           endif else output_cube(*,*,p)=D_OPLANE
        endfor 
     end
     (instfct eq 0): begin
        close, 2
        dlamb=dabs(lw-fw)/(nbw-1)*2.
        nbwref=nbw
        tab_wref = fltarr(nbwref)
        tab_dlamb = fltarr(nbwref)
        tab_index = fltarr(nbwref)

        for nbcw=0, nbw-1 do begin
           case dim1 of
              'mi': begin
                 case dim2 of
                    'mi': begin
                       tab_wref(nbcw)=tab_w(nbcw)
                       tab_dlamb(nbcw)=dlamb
                    end
                    'cm': begin
                                ;  change the unity of the spectral quantity used
                       tab_dlamb(nbw-nbcw+1)=1./tab_w(nbcw)^2*dlamb*1.d04
                       wave=1./tab_w(nbcw)*1.d04
                       tab_wref(nbw-nbcw+1)=wave
                    end
                    'nm': begin
                                ;  change the unity of the spectral quantity used
                       tab_wref(nbcw)=tab_w(nbcw)*1.d03
                       tab_dlamb(nbcw)=dlamb*1.d03
                    end
                    else: goto, mistake9105
                 endcase
              end
              'cm': begin
                 case dim2 of
                    'mi': begin
                                ;  change the unity of the spectral quantity used
                       tab_dlamb(nbw-nbcw+1)=1./tab_w(nbcw)^2*dlamb*1.d04
                       wave=1./tab_w(nbcw)*1.d04
                       tab_wref(nbw-nbcw+1)=wave
                    end
                    'cm' : begin
                       tab_wref(nbcw)=tab_w(nbcw)
                       tab_dlamb(nbcw)=dlamb
                    end
                    'nm': begin
                                ;  change the unity of the spectral quantity used
                       tab_dlamb(nbw-nbcw+1)=1./tab_w(nbcw)^2*dlamb*1.d07
                       wave=1./tab_w(nbcw)*1.d07
                       tab_wref(nbw-nbcw+1)=wave
                    end
                    else : goto, mistake9105
                 endcase
              end
              'nm' : begin
                 case dim2 of
                    'mi': begin
                       tab_wref(nbcw)=tab_w(nbcw)/1.d03 
                       tab_dlamb(nbcw)=dlamb/1.d03
                    end
                    'cm' : begin
                                ; change the unity of the spectral quantity used
                       tab_dlamb(nbw-nbcw+1)=1./tab_w(nbcw)^2*dlamb*1.d07
                       wave=1./tab_w(nbcw)*1.d07
                       tab_wref(nbw-nbcw+1)=wave
                    end
                    'nm': begin
                       tab_wref(nbcw)=tab_w(nbcw)
                       tab_dlamb(nbcw)=dlamb
                    end
                    else : goto, mistake9105
                 endcase
              end
              else : goto, mistake9100
           endcase
        endfor
        output_cube = fltarr(nbwref,qmax,pmax)
                                ;***Deresolution (actually just a copy with a change of the unit)
        for p=0, pmax-1 do begin
           for q=0, qmax-1 do begin
              readu, 3, refhp
              output_cube(*,q,p) = refhp
           endfor
        endfor 
        close, 3   
     end

  ;  (instfct eq 1): begin
    (instfct eq 4): begin  ;ATTENTION MODIF F. ANDRIEU
; CRISM case
        type='SWE'
        thermal_shift=0.
        case 1B of

           type eq 'SWE': begin
;typical CRISM parameters
              nbCol=604
              nbwref=437

;loading CRISM PSF's...
              print, 'loading CRISM PSFs ...'
              PSF_inst=fltarr(float(nbCol)*float(nbwref),10)
              PSF_SS=fltarr(nbwref,10)
              OpenR, lun, '/home/fandrieu/Inversion/Programmes/PSF_SS.txt', /Get_Lun
              ReadF, lun, PSF_SS
              Free_Lun, lun
              for j=0, nbwref-1 do for i=0, nbCol-1 do PSF_inst(j*float(nbCol)+i,*)=PSF_SS(j,*)
              PSF_inst=reform(PSF_inst,nbCol ,nbwref, 10, /overwrite)
; applying thermal shift to the entire PSF
              PSF_inst(*,*,3)=PSF_inst(*,*,3)+thermal_shift
              PSF_inst(*,*,6)=PSF_inst(*,*,6)+thermal_shift
              PSF_inst(*,*,9)=PSF_inst(*,*,9)+thermal_shift

;; dans le cas où la version desmilée de l'acquisition centrale est utilisée pour générer le CSP par binning
;;            nbCol=604
;;            nbwref=437
;;            PSF_inst_SS=fltarr(float(nbCol)*float(nbwref),10)
;;            PSF_SS=fltarr(nbwref,10)
;;            OpenR, lun, '/user/workdir/doute/Missions_inst/MRO/CRISM/calib/PSF_SS.txt', /Get_Lun
;;            ReadF, lun, PSF_SS
;;            Free_Lun, lun
;;            for j=0, nbwref-1 do for i=0, nbCol-1 do PSF_inst_SS(j*float(nbCol)+i,*)=PSF_SS(j,*)
;;            PSF_inst_SS=reform(PSF_inst_SS,nbCol ,nbwref, 10, /overwrite)
;; applying thermal shift to the entire PSF
;;            PSF_inst_SS(*,*,3)=PSF_inst(*,*,3)_SS+thermal_shift
;;            PSF_inst_SS(*,*,6)=PSF_inst(*,*,6)_SS+thermal_shift
;;            PSF_inst_SS(*,*,9)=PSF_inst(*,*,9)_SS+thermal_shift
              tab_wref=PSF_SS(*,6)/1.e3
              tab_dlamb=2.3548*sqrt(0.5/PSF_SS(*,5))/1.e3
           end

           (type eq 'FRT') or (type eq 'HRL') or (type eq 'CSP') : begin
;typical CRISM parameters
              nbCol=604
              nbwref=437

;loading CRISM PSF's...
              print, 'loading CRISM PSFs ...'
              PSF_inst=fltarr(float(nbCol)*float(nbwref),10)
              OpenR, lun, '/user/workdir/doute/Missions_inst/MRO/CRISM/calib/PSF_CRISM.txt', /Get_Lun
              ReadF, lun, PSF_inst
              Free_Lun, lun
              PSF_inst=reform(PSF_inst,nbCol ,nbwref , 10, /overwrite)
; applying thermal shift to the entire PSF
              PSF_inst(*,*,3)=PSF_inst(*,*,3)+thermal_shift
              PSF_inst(*,*,6)=PSF_inst(*,*,6)+thermal_shift
              PSF_inst(*,*,9)=PSF_inst(*,*,9)+thermal_shift

; applying to the channels that encompass the 2 micron CO2 band the experimental shift evaluated by the calibration procedure
              ;; for n_wref=140,172 do begin
              ;;    PSF_inst(2:601,n_wref,3)=PSF_inst(2:601,n_wref,3)+WA_flight_delta
              ;;    PSF_inst(2:601,n_wref,6)=PSF_inst(2:601,n_wref,6)+WA_flight_delta
              ;;    PSF_inst(2:601,n_wref,9)=PSF_inst(2:601,n_wref,9)+WA_flight_delta           
              ;; endfor
           end
        endcase
        tab_pt=ptrarr(nbwref) & tab_ind=lonarr(nbwref,2)
        output_cube = fltarr(nbwref,qmax,pmax)
        D_OPLANE=fltarr(nbwref,qmax)
        for p=0,pmax-1 do begin
           fct_deresol_block2,p,D_OPLANE,err_flag,err_msg,kresol=4, n_s1=218, n_s2=218
           if(err_flag ne 0) then begin
              goto, mistake9170
           endif else begin
              output_cube(*,*,p)=D_OPLANE
              print,format='("line ",I5," treated")',p
           endelse 
        endfor 
     end 

  endcase

  WAVE={BAND_BIN_ORIGINAL_BAND:indgen(nbwref), BAND_BIN_UNIT:dim2, BAND_BIN_CENTER:tab_wref, BAND_BIN_WIDTH:tab_dlamb}

  SOFTWARE={SOFTWARE_NAME:'SPECTRIBASE', SOFTWARE_DESC:'Produces synthetic spectral images cubes', SOFTWARE_VERSION_ID:'2012-01-16', SOFTWARE_RELEASE_DATE:'2012-01-16'}

  close, /all
                                ;create_cube, output_cube, 'BIP', final_output, CORE_NAME='NORMALIZED REFLECTANCE', WAVE=WAVE, SOFTWARE=SOFTWARE 
  create_cube, output_cube, 'BIP', final_output, CORE_NAME='NORMALIZED REFLECTANCE', WAVE=WAVE,DICO='/home/fandrieu/Simu_IPAG/pro/LPG_PDS_LABEL.dic', LABEL=label

  WRITE_ENVI_HEADER2, label, final_output, BNAMES = wave
  
  return
  
  mistake9000: begin
     print, 'Non existent environment file : ', 'env.txt'
     goto, mistake9999
  end
  
  mistake9095: begin
     print, 'Non existent spectral characteristic file : ', charaspectr
     goto, mistake9999
  end

  mistake9096: begin
     print, 'Non existent raw binary scratch file : ', scratch_file
     goto,  mistake9999
  end

  mistake9100: begin
     print, 'Invalid unit for the input optical files'
     goto, mistake9999
  end

  mistake9105: begin
     print, 'Invalid unit for the output cube file'
     goto, mistake9999
  end

  mistake9111: begin
     print, 'Invalid image option : ', opt3
     goto, mistake9999
  end

  mistake9112: begin
     print, 'Invalid correction option : ', opt4
     goto, mistake9999
  end

  mistake9113: begin
     print, 'Invalid interpolation option : ', opt5
     goto, mistake9999
  end

  mistake9150: begin
     print, 'Non existent deresolution file : ', derefich
     goto, mistake9999
  end

  mistake9170: begin
     print, 'Deresolution failed : ', err_msg
     goto, mistake9999
  end

  mistake9999: begin
  end


  
 ; stop

END


