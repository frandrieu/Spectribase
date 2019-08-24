  program spectribase

  use def_op_lst
  use conv_trigo

  implicit none

  !+
  !   $Id: spectribase.f90,v 1.2 2004/10/05 09:23:56 madec Exp doute $
  !
  ! NAME: 
  !	Spectribase (DOUSER subroutine under ISIS
  ! 	developing environment)
  !
  !
  !
  ! PURPOSE:
  !	 produces a synthetic spectral image cube 
  ! 	of a planetary surface given the chemical and physical 
  ! 	properties of the latter and the observation geometry 
  !
  !
  !
  ! CATEGORY: 
  !	hyperspectral imaging
  !
  !
  !
  ! CALLING SEQUENCE: 
  !	see Spectrimag manual
  !
  !
  !
  ! INPUTS: 
  !	see TAE and other parameters below  
  !
  !
  !
  ! OPTIONAL INPUTS:
  !
  !
  !
  ! KEYWORD PARAMETERS:
  !
  !
  !
  ! OUTPUTS: 
  !	an ISIS cube coding the hyperspectral image of the
  ! 	planetary surface. The physical organization of the cube is BSQ.
  ! 	Spatial dimensions of the image are given by the user. Spectral
  ! 	dimension is automatically calculated depending on the spectral
  ! 	characteristics of the instrument. Five backplanes are provided
  ! 	with the following information for each pixel of the image:
  ! 	mean latitude, longitude, incidence angle, emission angle 
  ! 	and phase angle.
  !
  !
  !
  ! OPTIONAL OUTPUTS:
  !
  !
  !
  ! COMMON BLOCKS:
  !    fctphase
  !
  !
  ! CALLS:
  !
  ! SIDE EFFECTS:
  !
  !
  !
  ! RESTRICTIONS:
  !
  !
  !
  ! PROCEDURE:
  !
  !
  !
  ! EXAMPLE:
  !
  !
  !
  ! MODIFICATION HISTORY:
  !	$Log: spectribase.f90,v $
  !	Revision 1.2  2004/10/05 09:23:56  madec
  !	*** empty log message ***
  !
  !	Revision 1.1  2004/10/01 09:10:53  madec
  !	Initial revision
  !   
  !
  !	27/02/2002 - Sylvain Doute' LPG - Version 7
  !	26/08/2013 - Francois Andrieu IDES - Photometric corrections
  !     28/11/2016 - François Andrieu GEOPS - Version Spectribase_2.1
  !-

  !+---------------------------------------------------------------------+
  !
  !                 ** DISORT I/O specifications **

  INTEGER, parameter::  MAXCLY = 20, MAXMOM = 299, MAXPHI = 180, MAXULV = 15, MAXUMU = 100
  character(len=127), parameter::  HEADER=''
  LOGICAL::  LAMBER, PLANK, ONLYFL, PRNT(5), USRANG, USRTAU
  INTEGER::  IBCND, NMOM, NLYR, NUMU, NUMU0, NSTR, NPHI, NTAU
  REAL::     ACCUR, ALBEDO, BTEMP, DTAUC( MAXCLY ), FBEAM, FISOT,&
       &TPHI( MAXPHI ), PMOM( 0:MAXMOM, MAXCLY ),&
       &PHI0, SSALB( MAXCLY ), TEMPER( 0:MAXCLY ), TEMIS, TTEMP,&
       &WVNMLO, WVNMHI, TUMU( MAXUMU ), TUMU0( MAXUMU ), UTAU( MAXULV ),&
       &TUMUS( MAXUMU ), TLEVEL( MAXULV), tab_alt( MAXULV) 

  REAL:: RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),&
       &DFDT( MAXULV ), UAVG( MAXULV ),&
       &UU( MAXUMU, MAXULV, MAXPHI ), ALBMED( MAXUMU ),&
       &TRNMED( MAXUMU )

  INTEGER:: rang(MAXUMU), tab_rang(1), ind_lyr( MAXULV)
  LOGICAL:: mask(MAXUMU)

  !+---------------------------------------------------------------------+

  !
  !------------------------------------------------------------

  !-------------------------------------------------------------
  !| TAE parameters                                            |
  !|         envfich file giving the relative paths of the     |
  !|         input text and binary files                       |
  !|         listarea file describing the area                 |
  !|           contributions                                   |
  !|         scenarea file describing the surface              |
  !|           representation of each area                     |
  !|         charactarea file describing the chemical and      |
  !|           physical characteristics of each area           |
  !|         geoarea file describing the geographical          |
  !|           extension of each area                          |
  !|         geopix file describing the geographical           |
  !|           distribution of the pixels                      |
  !|         geovisil file describing the local conditions     |
  !|           of illumination and observation                 |
  !|         charaspectr file describing the spectral          |
  !|           characteristics of the pixels and the way to    |
  !|           treat them                                      |
  !------------------------------------------------------------|

  !|-----------------------------------------------------------|
  !|  Other parameters :                                       |
  !|         name files of the compound optical indexes        |
  !|         {*dir} relative paths of the input files          |
  !|         homedir absolute path of the owner home directory |
  !|         opt1 area option                                  |
  !|         opt2 opposition effect option                     |
  !|         opt3 image option                                 |
  !|         opt4 correction option                            |
  !|         opt5 interpolation option                         |
  !|         opt6 optical depth option                         |
  !|         opt7 atmospheric correction option                |
  !|         contrib contribution flag of the current area     |
  !|         nbc nb of diff. components types                  |
  !|         name name of the file of the current              |
  !|           component                                       |
  !|         Bo amplitude of the opposition effect             |
  !|         wid angular width of the opposition effect        |
  !|         dzeta mean slope of the surface roughness (deg)   |
  !|         propsurf proportion of projected surface for the  |
  !|           current area (=1 for uniform area)              |
  !|         nbl nb of layers in the current area              |
  !|         nbatm nb of atmospheric layers for all area       |
  !|         tex textures of the diff. layers                  |
  !|         nbpc numbers of the present components            |
  !|           for the diff. layers                            |
  !|         nat orders of the components in each layers       |
  !|         tabcompa compactness of the layer matrix          |
  !|           in the case of compact textures                 |
  !|         tabthick real metric thicknesses of the layers    |
  !|           except the first one (infinite thick.)          |
  !|         tabwidth real metric thicknesses of the           |
  !|           atmospheric layers                              |
  !|           except the first one (infinite thick.)          |
  !|         tabg phase function parameters of the layers      |
  !|         tabprop proportions of the present                |
  !|           components in each layer                        |
  !|         tabdiam grain or inclusion sizes for each         |
  !|           layer (microns)                                 |
  !|         tabs internal scattering parameters for each      |
  !|           layer (microns-1)                               |
  !|         dim1 unit of the spectral quantity used in        |
  !|           optical files                                   |
  !|         fww first wave of the working domain within the   |
  !|            optical files                                  |
  !|         lww last wave of the working domain within the    |
  !|            optical files                                  |
  !|         dim2 unit of the spectral quantity used for the   |
  !|           spectral image cube                             |
  !|         pmax vertical dimension of the spectral           |
  !|           image cube                                      |
  !|         qmax horizontal dimension of the spectral         |
  !|            image cube                                     |
  !|         instfct type of instrument spectral response      |
  !|           (gaussian, triangular, trapezoidal)             |
  !|         dlamb instrument spectral response FWHM           |
  !|         topdlamb second FWHM for a trapezoidal spectral   |
  !|            response function                              |
  !|         slamb final spectral sampling after the           |
  !|           deresolution                                    |
  !|         derefich file giving the reference wavelengths    |
  !|           from which the deresolution is made             |
  !|         derefichdir relative path of the derefich         |
  !|           directory                                       |
  !|         firstw first wave of the convolved spectra        |
  !|         lastw last wave of the convolved spectra          |
  !|         u latitude rank in the integration grid-point     |
  !|           system                                          |
  !|         v longitude rank in the integration grid-point    |
  !|           system                                          |
  !|         q pixel horizontal rank in the spect. imag. cube  |
  !|         p pixel vertical rank in the spect. imag. cube    |
  !|         latmin lower latitude (deg) of the geographical   |
  !|           integration grid                                |
  !|         latmax upper latitude (deg) of the geographical   |
  !|           integration grid                                |
  !|         umax nb of grid points along the longitude        |
  !|           lines                                           |
  !|         longmin lower longitude (deg) of the geographical |
  !|           integration grid                                |
  !|         longmax upper longitude (deg) of the geographical |
  !|           integration grid                                |
  !|         vmax nb of grid points along the reference        |
  !|           latitude line                                   |
  !|         muo cosine of the incident zenithal angle         |
  !|            corrected if surface roughness taken into      |
  !|            account                                        |
  !|         mu cosine of the emergent zenithal angle          |
  !|            corrected if surface roughness taken into      |
  !|            account                                        |
  !|         phi azimuthal emergent angle (rad)                |
  !|         alpha local scattering angle (deg)                | 
  !|         wave wavelength or wavenumber (microns, nanometers|
  !|          or cm-1)                                         |
  !|         tabindreal real refractive indexes of all the     |
  !|           component types                                 |
  !|         tabindimag imaginary refractive indexes           |
  !|           of all the component types                      |
  !|         fctP real phase function                          |
  !|         g real phase function scattering parameter        |
  !|         nbcmax maximum number of components               |
  !|           per layer                                       |
  !|         nbord number of Legendre components taking        |
  !|           into account in radiative transfer              |
  !|           calculations                                    |
  !|         nbordm maximum value of the previous quantity     |
  !|         eps infinitesimal quantity for the calculation    |
  !|           of the limits                                   |
  !|                                                           |
  !|-----------------------------------------------------------|

  !|-----------------------------------------------------------|
  !|  Partial results :                                        |
  !|         direct access scratch file linear array           |
  !|           containing the sequential list of all the       |
  !|           spectro-pixels                                  |
  !|         pix spectral contribution of the current area     |
  !|           to the pixels                                   |
  !|         pixa contribution to the pixels of a              |
  !|           lambertian surface corresponding to the         |
  !|           current area                                    |
  !|         refl spectral contribution of the current         |
  !|           geographical point to the pixel involved        |
  !|         ap matrices of Legendre coefficients              |
  !|         tabomeg single scattering albedos of the          |
  !|            layers of the current area                     |
  !|         vm nb of grid points along the latitude           |
  !|           line u                                          |
  !|         nbw number of wavelengths or wavenumbers for the  |
  !|           raw spectra                                     |
  !|         These wavelengths or wavenumbers are stored       |
  !|           within the allocatable array tab_w              | 
  !|         nbwref number of wavelengths or wavenumbers       |
  !|           for the convolved spectra                       |
  !|                                                           |
  !|-----------------------------------------------------------|

  !|-----------------------------------------------------------|
  !|  Functions :                                              |
  !|         singalbg single scattering albedo for a           |
  !|           granular texture                                |
  !|         singalbc single scattering albedo for a           |
  !|           compact texture                                 |
  !|         fctP real phase function                          |
  !|         derefct deresolution function                     | 
  !|         rad transformation deg->rad                       |
  !|                                                           |
  !|-----------------------------------------------------------|

  !|-----------------------------------------------------------|
  !|                                                           |
  !|  Other quantities :                                       |
  !|         up rank of the previous latitude line in the      |
  !|           area description                                |
  !|         vp rank of the previous longitude point in        |
  !|           the area description                            |
  !|         nbcw rank of the current wavelength or            |
  !|           wavenumber                                      |
  !|         nb index of component types                       |
  !|         nbta number of current treated area               |
  !|         k index of present components                     |
  !|         m index of layers                                 |
  !|         n index of Legendre order                         |
  !|         i,j,t loop indexes                                |
  !|         dlat latitude step                                |
  !|         dlong longitude step                              |
  !|         angl any angle                                    |
  !|         omeg current single scattering albedo             |
  !|         tau current optical depth                         |
  !|         surf surface of the integration element           |
  !|         surfref surface of the integration element        |
  !|           at the reference latitude                       |
  !|         thetauo incident zenithal angle                   |
  !|         thetau emergent zenithal angle                    |
  !|         fw first wave of the index files                  |
  !|         lw last wave of the index files                   |
  !|         unit1-2 indicators of spectral units              |
  !|         str_unit2 string containing the unit of the       |
  !|            output spectral cube wavelengths               |
  !|         str_length2 its corresponding length              |
  !|         *ptm addresses of chain list heads                |
  !|         eta conversion variable                           |
  !|                                                           |
  !|-----------------------------------------------------------|

  !-----------------------------------------------------------------
  !  Definition of internal types and variables

  integer(kind=2):: u,v,umax,vmax,vm,up,vp,q
  integer(kind=2), save:: p,pmax,qmax
  integer(kind=2):: i,j,t,m,k,n,nbc,nb,r
  integer(kind=2):: nbta,nb_digit
  integer(kind=2), save :: unit1, unit2, instfct
  integer(kind=2):: nbc_block, nb_sfiles, sf_nb
  integer(kind=2), save:: nb_block
  integer(kind=2):: pos1,pos2, bool_rough  

  integer(kind=2), parameter:: nbcmax=15

  integer(kind=4):: nb_args, status, ERRVAL
  integer(kind=4):: nbl,nbord, nbatm
  integer(kind=4):: nbcw,nbcwp, rec_nb,count_nb_k, count_nb_area
  integer(kind=4), save :: nbw,nbwref
  integer(kind=4):: err_flag, err_stat, ios, selector
  integer(kind=4):: I_RECLENGTH
  integer(kind=4):: str_length2
  integer(kind=4):: nb_pt_skip,n_pt_skip
  integer(kind=4),parameter:: nb_pt_skip_max=25
  integer(kind=4), save:: nbwp, recpmax
  integer,parameter :: spherical=0 !1 to compute spherical lambertian reflectance
! ; anything else to compute Hapke model in case of granular monolayer
  integer(kind=2), dimension(:), allocatable :: nbpc
  integer(kind=4), dimension(:), allocatable :: expo
  integer(kind=2), dimension(:,:), allocatable :: nat

  real(kind=4):: latmin,latmax,latm,dlat,dlong,longmin,longmax
  real(kind=4):: muo,mu,thetauo,thetau,phi,alpha,Bo,wid
  real(kind=4):: muo_a,mu_a,muo_b,mu_b,phi_a
  real(kind=4):: depth, g, theta, theta0, fct_reflectance
  real(kind=4):: betap, betag, fact_D, altitude, sumtau
  real(kind=4):: indrs, indis, R0p, T0p, R_diff, R_spec
  real(kind=4):: gamma, ru, rl, rc, tc, rs, T_diff
  real(kind=4):: fact, sep, coefreflgp, est_fact
  real(kind=4):: fct_reflectance_SC_spec, fct_reflectance_SC_diff
  real(kind=4),save :: dlamb,topdlamb
  real(kind=4) :: slamb,firstw,lastw,fct_S,dzeta, dzeta_old
  real(kind=4), parameter:: pi=3.141592653589793
  real(kind=4):: delta_R_moy
  real(kind=4), parameter :: delta_R_max=1.0e-2 
  real(kind=4), dimension(:,:), pointer, save:: pix_mem
  real(kind=4), dimension(:,:,:), allocatable, save:: pix_bck
  real(kind=4), dimension(:), pointer,save:: tab_w
  real(kind=4), dimension(:,:), allocatable:: tauatm

  real(kind=8):: wave,omeg,b,pixm,tau,propsurf,singalbg,singalbc, dn, dk, dn2, dk2
  real(kind=8):: refl, R0, T0, se, si,last_tabindreal,last_tabindimag, last_indre, last_indim
  real(kind=8):: fww,lww,fw,lw, perc_wave, perc_mem, perc_comp


  real(kind=8), parameter:: eps=1.d-4
  real(kind=4), parameter:: mem_max1=1.d9 !previously 5.d6 ! size limit in memory (in octet)
  real(kind=4), parameter:: mem_max2=2.d9 ! file size limit on system (in octet) : maximum size of scratch file = 2Go

  real(kind=8),dimension(:),allocatable:: tabcompa,tabtau,tabthick,tabg, tabwidth
  real(kind=8),dimension(:,:),allocatable:: tabprop,tabdiam,tabs
  real(kind=8),dimension(:),allocatable:: tabindreal,tabindimag,tabdindreal,tabdindimag, tabn, tabk
  real(kind=4),dimension(:,:),allocatable:: ap
  real(kind=8),dimension(:),pointer:: cfa,cfb
  real(kind=8),dimension(:),allocatable:: indr,indi,prop,diam
  real(kind=8),dimension(:),allocatable:: s,tabomeg
  real(kind=8),dimension(:,:),pointer,save:: pix,pixa
  real(kind=8), dimension(:,:), allocatable:: tab_fact, tab_sep, tab_se
  real(kind=8), dimension(:,:), allocatable:: tab_spec

  character(len=255):: err_msg
  character(len=255):: listarea,scenarea,charactarea,raw_output
  character(len=255):: geoarea,geopix,geovisil,atmtau,name
  character(len=255):: charaspectr,fichspect,derefich, scratch_file
  character(len=200):: suffix
  character(len=90):: namedir
  character(len=90):: homedir
  character(len=90):: envfich,listareadir,scenareadir
  character(len=90):: charactareadir,geoareadir,geopixdir,geovisildir
  character(len=90):: charaspectrdir,derefichdir, atmtaudir
  character(len=20):: fmt_string, ERRKEY
  character(len=10):: str_unit2, name_prog
  character(len=4):: opt1,opt3
  character(len=2):: dim1,dim2
  character(len=1):: contrib,opt2,opt4,opt5,opt6,opt7
  character(len=2):: err_code_mem

  character(len=1),dimension(:),allocatable:: tex

  LOGICAL, parameter:: back=.true.

  !  type(DXML_D_FCT_STRUCTURE):: fct_struct

  type(wnode), pointer:: w_item
  type(indnode), pointer:: ind_item
  type(geo1node), pointer:: geo1_item
  type(geo2node), pointer:: geo2_item
  type(wpt):: wptm 
  type(indpt), dimension(:), pointer:: indptm_tab
  type(indpt), dimension(:), pointer:: indptm
  type(geo1pt):: geo1ptm 
  type(geo2pt), dimension(:), pointer:: geo2ptm_tab

  common/fctphase/omeg,b,tau

  !+---------------------------------------------------------------------+

  !     .. External Subroutines ..

  EXTERNAL  DISORT, GETMOM, ERRMSG

  !     ..
  !     .. Intrinsic Functions ..

  INTRINSIC ASIN, FLOAT, INDEX

  !
  !-----------------------------------------------------------------

  !*** Read the input file names




  nb_args=NARGS()
  if(nb_args.ne.12) then
     write(unit=6,FMT='(A255)')'Syntax : spectribase TO="out.cub"  SFROM=" " ENVFICH="env.txt" LISTAREA="list_suffix.txt" SCENAREA="scen_suffix.txt" CHARACTA="chara_suffix.txt" GEOAREA="geoa_suffix.bin" GEOPIX="geop_suffix.bin" GEOVISIL="geov_suffix.bin" CHARASPE="spect_suffix_C.txt" USERNOTE=" "'
     stop
  endif

  write(err_msg,FMT='("******>> Starting the simulation <<******")')
  WRITE(UNIT=6,FMT='(A255)')err_msg

  write(err_msg,FMT='("**** Reading basic information *****")')
  WRITE(UNIT=6,FMT='(A255)')err_msg

  call getarg(0,name_prog,status)

  call getarg(1, raw_output, status)
  if (status.eq.-1) goto 9010
  raw_output=raw_output(4:len(raw_output)-1)
  write(unit=6,FMT='(1X,"Output file : ",A255)')raw_output

  CALL getarg(3, envfich, status)
  if (status.eq.-1) goto 9010
  envfich=envfich(9:len(envfich)-1)

  CALL getarg(4, listarea, status)
  listarea=listarea(10:len(listarea)-1)
  if (status.eq.-1) goto 9015

  CALL getarg(5, scenarea, status)
  scenarea=scenarea(10:len(scenarea)-1)
  if (status.eq.-1) goto 9020

  CALL getarg(6, charactarea, status)
  charactarea=charactarea(10:len(charactarea)-1)
  if (status.eq.-1) goto 9025

  CALL getarg(7, geoarea, status)
  geoarea=geoarea(9:len(geoarea)-1)
  if (status.eq.-1) goto 9030

  CALL getarg(8, geopix, status)
  geopix=geopix(8:len(geopix)-1)
  if (status.eq.-1) goto 9035

  CALL getarg(9, geovisil, status)
  geovisil=geovisil(10:len(geovisil)-1)
  if (status.eq.-1) goto 9040

  CALL getarg(10, charaspectr, status)
  charaspectr=charaspectr(10:len(charaspectr)-1)
  if (status.eq.-1) goto 9045

  !***  building the suffix array of characters
  pos1=index(charaspectr, 'spect')
  pos2=index(charaspectr, '.txt', back)
  suffix=charaspectr(pos1+5:pos2-1)
  write(err_msg,FMT='("Suffix of the experiment : ",A200)')suffix(2:)
  WRITE(UNIT=6,FMT='(A255)')err_msg

  !*** TO DO optional argument atmtau to be implemented if needed
  atmtau='none'
  if ((status.eq.-1).or.(trim(adjustl(atmtau)).eq.'none')) then 
     opt7='n' 
  else 
     opt7='y'
  endif

  !*** Read the input file relative paths

  open(unit=10,file=envfich,status="old",iostat=ios)
  if(ios/=0) goto 9000

  read(unit=10,FMT='(28X,A90)',iostat=ios) namedir

  read(unit=10,FMT='(17x,A90)',iostat=ios) homedir

  read(unit=10,FMT='(21x,A90)',iostat=ios) listareadir

  read(unit=10,FMT='(21X,A90)',iostat=ios) scenareadir

  read(unit=10,FMT='(24X,A90)',iostat=ios) charactareadir

  read(unit=10,FMT='(20X,A90)',iostat=ios) geoareadir

  read(unit=10,FMT='(19X,A90)',iostat=ios) geopixdir

  read(unit=10,FMT='(21X,A90)',iostat=ios) geovisildir

  read(unit=10,FMT='(21X,A90)',iostat=ios) atmtaudir

  read(unit=10,FMT='(24X,A90)',iostat=ios) charaspectrdir

  read(unit=10,FMT='(30x,A90)',iostat=ios) derefichdir
  if(ios/=0) goto 9005

  close(10)

  homedir=trim(adjustl(homedir))
  namedir=trim(adjustl(namedir))

  derefichdir=trim(adjustl(derefichdir))
  raw_output=trim(adjustl(raw_output))
  write(fmt_string,FMT='("(A",I2.2,",A",I2.2,",A1",",A",I2.2,")")'),&
       &len_trim(homedir),len_trim(derefichdir),len_trim(raw_output)
  write(raw_output,FMT=fmt_string)homedir,derefichdir,'/',raw_output

  listareadir=trim(adjustl(listareadir))
  listarea=trim(adjustl(listarea))
  write(fmt_string,FMT='("(A",I2.2,",A",I2.2,",A1",",A",I2.2,")")'),&
       &len_trim(homedir),len_trim(listareadir),len_trim(listarea)
  write(listarea,FMT=fmt_string)homedir,listareadir,'/',listarea

  scenareadir=trim(adjustl(scenareadir))
  scenarea=trim(adjustl(scenarea))
  write(fmt_string,FMT='("(A",I2.2,",A",I2.2,",A1",",A",I2.2,")")'),&
       &len_trim(homedir),len_trim(scenareadir),len_trim(scenarea)
  write(scenarea,FMT=fmt_string)homedir,scenareadir,'/',scenarea

  charactareadir=trim(adjustl(charactareadir))
  charactarea=trim(adjustl(charactarea))
  write(fmt_string,FMT='("(A",I2.2,",A",I2.2,",A1",",A",I2.2,")")'),&
       &len_trim(homedir),len_trim(charactareadir),len_trim(charactarea)
  write(charactarea,FMT=fmt_string)homedir,charactareadir,'/',charactarea

  geoareadir=trim(adjustl(geoareadir))
  geoarea=trim(adjustl(geoarea))
  write(fmt_string,FMT='("(A",I2.2,",A",I2.2,",A1",",A",I2.2,")")'),&
       &len_trim(homedir),len_trim(geoareadir),len_trim(geoarea)
  write(geoarea,FMT=fmt_string)homedir,geoareadir,'/',geoarea

  geopixdir=trim(adjustl(geopixdir))
  geopix=trim(adjustl(geopix))
  write(fmt_string,FMT='("(A",I2.2,",A",I2.2,",A1",",A",I2.2,")")'),&
       &len_trim(homedir),len_trim(geopixdir),len_trim(geopix)
  write(geopix,FMT=fmt_string)homedir,geopixdir,'/',geopix

  geovisildir=trim(adjustl(geovisildir))
  geovisil=trim(adjustl(geovisil))
  write(fmt_string,FMT='("(A",I2.2,",A",I2.2,",A1",",A",I2.2,")")'),&
       &len_trim(homedir),len_trim(geovisildir),len_trim(geovisil)
  write(geovisil,FMT=fmt_string)homedir,geovisildir,'/',geovisil

  atmtaudir=trim(adjustl(atmtaudir))
  atmtau=trim(adjustl(atmtau))
  write(fmt_string,FMT='("(A",I2.2,",A",I2.2,",A1",",A",I2.2,")")'),&
       &len_trim(homedir),len_trim(atmtaudir),len_trim(atmtau)
  write(atmtau,FMT=fmt_string)homedir,atmtaudir,'/',atmtau

  charaspectrdir=trim(adjustl(charaspectrdir))
  charaspectr=trim(adjustl(charaspectr))
  write(fmt_string,FMT='("(A",I2.2,",A",I2.2,",A1",",A",I2.2,")")'),&
       &len_trim(homedir),len_trim(charaspectrdir),len_trim(charaspectr)
  write(charaspectr,FMT=fmt_string)homedir,charaspectrdir,'/',charaspectr

  !* initialization of the number of current area treated 
  nbta=1

  !*************************************************************
  !*** Opening the image file

  open(14,file=charaspectr,status="old",iostat=ios)
  if(ios/=0) goto 9095

  read(14,FMT='(54x,A2)',iostat=ios)dim1
  if(ios/=0) goto 9110
  select case(trim(adjustl(dim1)))
  case('mi')
     unit1=2
  case('cm')
     unit1=1
  case('nm')
     unit1=3
  case default
     goto 9100
  end select

  if(dim1.eq.'mi') then
     read(14,FMT='(48x,f7.4,1x,f7.4)',iostat=ios)fww,lww
  else
     read(14,FMT='(48x,f7.1,1x,f7.1)',iostat=ios)fww,lww
  endif

  read(14,FMT='(52x,A2)',iostat=ios)dim2
  if(ios/=0) goto 9110
  select case(trim(adjustl(dim2)))
  case('mi')
     unit2=2
     str_unit2='MICROMETER'
     str_length2=10
  case('cm')
     unit2=1
     str_unit2='WAVENUMBER'
     str_length2=10
  case('nm')
     unit2=3
     str_unit2='NANOMETER'
     str_length2=9
  case default
     goto 9105
  end select

  read(14,FMT='(42x,I5)',iostat=ios)pmax
  write(err_msg,FMT='("Nb of image lines : ",I5)')pmax
  WRITE(UNIT=6,FMT='(A255)')err_msg

  read(14,FMT='(44x,I5)',iostat=ios)qmax
  write(err_msg,FMT='("Nb of image samples : ",I5)')qmax
  WRITE(UNIT=6,FMT='(A255)')err_msg

  if(ios/=0) goto 9110

  !*************************************************************

  !*** Initialization of the pixel tables

  allocate(pix(qmax,pmax),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='12' 
     goto 9065
  end if
  allocate(pixa(qmax,pmax),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='13' 
     goto 9065
  end if
  allocate(pix_bck(qmax,pmax,7),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='14b' 
     goto 9065
  end if

  !* initialization to zero of the pixel tables
  pixa=0.
  pix=0.
  pix_bck=0.

  !*************************************************************

  !*** Initialization of the general geographical dynamic structure

  open(10,file=geovisil,form='unformatted',status="old",iostat=ios)
  if(ios/=0) goto 9050

  open(11,file=geopix,form='unformatted',status="old",iostat=ios)
  if(ios/=0) goto 9055

  !* description of the work geographical space
  read(10,iostat=ios)latmin,latmax,umax,longmin,longmax,&
       &vmax
  if(ios/=0) goto 9060

  write(err_msg,FMT='("Geographical limits and grid points")')
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  latmin : ",f5.1)')latmin
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  latmax : ",f5.1)')latmax
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  umax : ",I5)'),umax
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  longmin : ",f6.1)')longmin
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  longmax : ",f6.1)')longmax
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  vmax : ",I4)')vmax
  WRITE(UNIT=6,FMT='(A255)')err_msg

  !* allocation of the array containing the pointers to the 
  !latitude chain lists
  allocate(geo2ptm_tab(umax),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='1' 
     goto 9065
  end if

  dlat=(latmax-latmin)/umax

  !* calculation of the reference latitude
  if(abs(latmax).gt.abs(latmin)) then
     latm=abs(latmax)
  else
     latm=abs(latmin)
  endif

  loop_lat: do u=1,umax

     dlong=(longmax-longmin)/vmax

     call init_geo2(geo2ptm_tab(u)%head,geo2ptm_tab(u)%tail)

     allocate(geo2_item, stat=err_stat)
     ! check if the allocation is successfull 
     if(err_stat/=0) then
        err_code_mem='2' 
        goto 9065
     end if

     loop_long: do v=1,vmax
        read(10,iostat=ios)geo2_item%muo,geo2_item%mu,&
             &geo2_item%phi,geo2_item%alpha, geo2_item%depth

        if(ios<0) exit loop_lat
        if(ios>0) goto 9060
        read(11,iostat=ios)geo2_item%p,geo2_item%q
        if(ios<0) exit loop_lat
        if(ios>0) goto 9070

        ! Calculation of the mean conditions of illumination and viewing for each pixel
        pix_bck(geo2_item%q,geo2_item%p,1)=pix_bck(geo2_item%q,geo2_item%p,1)&
             &+latmin+(FLOATI(u)-0.5)*dlat
        pix_bck(geo2_item%q,geo2_item%p,2)=pix_bck(geo2_item%q,geo2_item%p,2)&
             &+longmin+(FLOATI(v)-0.5)*dlong
        pix_bck(geo2_item%q,geo2_item%p,3)=pix_bck(geo2_item%q,geo2_item%p,3)&
             &+acos(geo2_item%muo)/pi*180.
        pix_bck(geo2_item%q,geo2_item%p,4)=pix_bck(geo2_item%q,geo2_item%p,4)&
             &+acos(geo2_item%mu)/pi*180.
        pix_bck(geo2_item%q,geo2_item%p,5)=pix_bck(geo2_item%q,geo2_item%p,5)&
             &+(180.-geo2_item%alpha)
        pix_bck(geo2_item%q,geo2_item%p,6)=pix_bck(geo2_item%q,geo2_item%p,6)&
             &+geo2_item%depth
        pix_bck(geo2_item%q,geo2_item%p,7)=pix_bck(geo2_item%q,geo2_item%p,7)+1.


        call add_geo2(geo2_item,geo2ptm_tab(u)%head,geo2ptm_tab(u)%tail)
        allocate(geo2_item, stat=err_stat)
        ! check if the allocation is successfull 
        if(err_stat/=0) then
           err_code_mem='3' 
           goto 9065
        end if

     end do loop_long

  end do loop_lat

  write(err_msg,FMT='("Reading files geovisil and geopix : OK")')
  WRITE(UNIT=6,FMT='(A255)')err_msg

  close(10)
  close(11)

  !*************************************************************

  !*** Opening the area files

  open(10,file=scenarea,status="old",iostat=ios)
  if(ios/=0) goto 9075
  open(11,file=charactarea,status="old",iostat=ios)
  if(ios/=0) goto 9080
  open(12,file=geoarea,form='unformatted',status="old",iostat=ios)
  if(ios/=0) goto 9085
  open(13,file=listarea,status="old",iostat=ios)
  if(ios/=0) goto 9090

  !*************************************************************

  !*** Initialization of the refractive index chain lists 

  read(10,FMT='(34x,I2)',iostat=ios)nbc
  if(ios/=0) goto 9115

  allocate(indptm_tab(nbc),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='4' 
     goto 9065
  end if
  allocate(indptm(nbc),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='5' 
     goto 9065
  end if
  allocate(tabindreal(nbc),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='6' 
     goto 9065
  end if
  allocate(tabindimag(nbc),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='7' 
     goto 9065
  end if

  nbw=0

  !* initialization of the wavelength chain list
  call init_w(wptm%head,wptm%tail)
  allocate(w_item,stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='8' 
     goto 9065
  end if

  loop_comp: do nb=1,nbc

     read(10,FMT='(3x,A20)',iostat=ios)name
     if(ios/=0) goto 9115

     name=trim(adjustl(name))
     write(fmt_string,FMT='("(A",I2.2,",A1,A",I2.2,")")'),&
          &len_trim(namedir),len_trim(name)
     write(name,FMT=fmt_string)namedir,'/',name
     open(15,file=name,form='formatted',status="old",iostat=ios)
     if(ios/=0) goto 9120

     !* filling the index chain list for the current compound
     call init_ind(indptm_tab(nb)%head,indptm_tab(nb)%tail)
     allocate(ind_item,stat=err_stat)
     ! check if the allocation is successfull
     if(err_stat/=0) then
        err_code_mem='9' 
        goto 9065
     end if

     do

        read(15,*,iostat=ios)w_item%wave,ind_item%indreal,&
             &ind_item%indimag
        if(ios<0) exit
        if(ios>0) goto 9125

        if((w_item%wave.ge.fww).and.(w_item%wave.le.lww)) then
           if(ind_item%indreal<1) print*,'indice reel non valable',' ',name
           if(ind_item%indimag<0) print*,'indice imaginaire non valable',' ',name 
           !* counting the number of wavelengths or wave numbers
           if(nb.eq.1) then
              if(nbw.eq.0) fw=w_item%wave
              nbw=nbw+1
              call add_w(w_item,wptm%head,wptm%tail)
              allocate(w_item,stat=err_stat)
              ! check if the allocation is successfull
              if(err_stat/=0) then
                 err_code_mem='10' 
                 goto 9065
              end if
           endif
           call add_ind(ind_item,indptm_tab(nb)%head,indptm_tab(nb)%tail)
           allocate(ind_item,stat=err_stat)
           ! check if the allocation is successfull
           if(err_stat/=0) then
              err_code_mem='11' 
              goto 9065
           end if
        endif

     end do

     close(15)

  end do loop_comp

  write(err_msg,FMT='("Reading optical constant files : OK")')
  WRITE(UNIT=6,FMT='(A255)')err_msg


  !*************************************************************

  ! Opening and initialization of the direct access scratch file containing the
  ! spectral image cube being built. Physical storage is BIP and data are stored 
  ! by blocks of mem_max1.

  ! determination of the block number
  nb_block=0
  do
     if((nb_block*mem_max1)>=(nbw*4.*qmax)) exit
     nb_block=nb_block+1
  enddo
  nbwp=nbw/nb_block ! determination of the max. index of each block
  if(((real(nbw)/real(nb_block))-nbwp) > tiny(1.)) nb_block=nb_block+1
  ! printing the block number and max index
  write(err_msg,FMT='("Number of blocks used: ",I4)') nb_block
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("maximum index of each block: ",I6)') nbwp
  WRITE(UNIT=6,FMT='(A255)')err_msg     

  ! allocation of the block image in memory
  allocate(pix_mem(nbwp,qmax),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='14a' 
     goto 9065
  end if
  pix_mem=0.


  allocate(tabdindreal(nbwp-1),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='6' 
     goto 9065
  end if
  allocate(tabdindimag(nbwp-1),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='7' 
     goto 9065
  end if
  allocate(tabn(nbwp),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='6' 
     goto 9065
  end if
  allocate(tabk(nbwp),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='7' 
     goto 9065
  end if


  ! determination of the scratch file number
  nb_sfiles=0
  do
     if((nb_sfiles*mem_max2)>=(real(pmax)*real(qmax)*real(nb_block)*real(nbwp)*4.)) exit
     nb_sfiles=nb_sfiles+1
  enddo
  ! determination of the max. record index for each scratch file
  recpmax=pmax*nb_block/nb_sfiles
  if(((real(pmax*nb_block)/real(nb_sfiles))-recpmax) > tiny(1.)) nb_sfiles=nb_sfiles+1
  ! printing the block number and max index
  write(err_msg,FMT='("Disk space used for the scratch file(s): ",f8.2," Mo")')&
       &(real(pmax)*real(qmax)*real(nb_block)*real(nbwp)*4./1.e06)
  WRITE(UNIT=6,FMT='(A255)')err_msg     
  write(err_msg,FMT='("Number of scratch files needed: ",I4)')nb_sfiles
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("maximum record index for each scratch file: ",I6)')recpmax
  WRITE(UNIT=6,FMT='(A255)')err_msg

  I_RECLENGTH=qmax*nbwp

  do k=1,nb_sfiles
     if (k.lt.9)  write(fmt_string,FMT='("(A",I2.2,",A8,I1)")'),&
          &len_trim(suffix)
     if ((k.ge.10).and.(k.lt.99)) write(fmt_string,FMT='("(A",I2.2,",A8,I2)")'),&
          &len_trim(suffix)
     write(scratch_file,FMT=fmt_string)suffix,'_scratch',k
     open(19+k, file=scratch_file, defaultfile='/data/fschmidt/fandrieu/fastworkdir', status='NEW', &
          &access='DIRECT', RECL=I_RECLENGTH, iostat=ios) 
     if(ios/=0) goto 9175
  end do

!  do p=1,pmax
!     do nbc_block=1,nb_block
!        sf_nb=((p-1)*nb_block+nbc_block-1)/recpmax+1
!        rec_nb=mod((p-1)*nb_block+nbc_block-1,recpmax)+1
!        write(UNIT=19+sf_nb,REC=rec_nb,iostat=ios)pix_mem
!    end do
!  end do
! if(ios/=0) goto 9176

  scratch_file = trim(adjustl(suffix))//'_scratch_header'
  open(19, file=scratch_file, form='unformatted', defaultfile='/data/fschmidt/fandrieu/fastworkdir', status='NEW', iostat=ios)
  if(ios/=0) goto 9175
  write(UNIT=19,iostat=ios)nb_block,nbwp,nb_sfiles,recpmax,nbw 
  if(ios/=0) goto 9176

  !*************************************************************

  !*** Read the different options

  !* area option
  read(10,*,iostat=ios)
  read(10,FMT='(14x,A4)',iostat=ios)opt1
  opt1=trim(adjustl(opt1))
  if((opt1/='diff').and.(opt1/='iden')) goto 9116

  !* optical depth option
  read(10,FMT='(23x,A1)',iostat=ios)opt6
  opt6=trim(adjustl(opt6))
  if((opt6/='y').and.(opt6/='n')) goto 9119

  if(ios/=0) goto 9115

  ! * nb of atmospheric layers
  read(10,FMT='(23x,I2)',iostat=ios)nbatm

  if(ios/=0) goto 9115

  !* image option
  read(14,FMT='(15x,A4)',iostat=ios)opt3
  opt3=trim(adjustl(opt3))
  if((opt3/='raw').and.(opt3/='conv')) goto 9111

  !* correction option
  read(14,FMT='(20x,A1)',iostat=ios)opt4
  opt4=trim(adjustl(opt4))
  if((opt4/='n').and.(opt4/='s').and.(opt4/='d')) goto 9112

  !* interpolation option
  read(14,FMT='(23x,A1)',iostat=ios)opt5
  opt5=trim(adjustl(opt5))
  if((opt5/='y').and.(opt5/='n')) goto 9113

  !* number of Fourier components used
  read(14,FMT='(36x,I3)',iostat=ios)nbord
  allocate(expo(nbord+1),stat=err_stat)
  expo=(/ (i,i=0,nbord) /)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='56' 
     goto 9065
  end if

  if(ios/=0) goto 9110

  write(err_msg,FMT='("Options summary : ")')
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  area option : ",A4)')opt1
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  optical depth option : ",A1)')opt6
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  image option : ",A4)')opt3
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  correction option : ",A1)')opt4
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  interpolation option : ",A1)')opt5
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  number of Fourier comp. : ",I2)')nbord
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("  Atmospheric correction. : ",A1)')opt7
  WRITE(UNIT=6,FMT='(A255)')err_msg
  write(err_msg,FMT='("Treatment of the different area")')
  WRITE(UNIT=6,FMT='(A255)')err_msg

  write(err_msg,FMT='("**** Starting the calculation *****")')
  WRITE(UNIT=6,FMT='(A255)')err_msg


  !-----------------------------------------------------------------
  !  Configuration of the calls to the radiative transfer engine (DISORT)

  PRNT(1)=.FALSE.
  PRNT(2)=.FALSE.
  PRNT(3)=.FALSE.
  PRNT(4)=.FALSE.
  PRNT(5)=.FALSE.
  NSTR = 16 ! number of computational polar angles
  USRTAU    = .TRUE. ! radiant quantities are to be returned at 
  USRANG    = .TRUE.! radiant quantities are to be returned at user specified polar angles  
  IBCND     = 0 ! general case
  FISOT      = 0.0 ! isotropic illumination at 0
  LAMBER    = .TRUE. ! isotropically reflecting boundary layer
  ONLYFL    = .FALSE. ! return intensities
  ALBEDO    = 0.0 ! albedo of the reflecting boundary layer
  PLANK     = .FALSE. ! ignore all thermal emission
  ACCUR = 0.0 ! convergence criteria
  PHI0      = 0.0 ! azimuth angle of incident beam

  !
  !-----------------------------------------------------------------


  !*************************************************************
  !*************************************************************
  !**                                                         **
  !**  Preparation of the data related to the first area      **
  !**                                                         **
  !*************************************************************
  !*************************************************************

  !*** Initialization of the layer characteristic arrays
  !    for the first area

  read(10,*,iostat=ios)
  read(10,*,iostat=ios)
  read(10,*,iostat=ios)
  read(10,*,iostat=ios)
  if(ios/=0) goto 9115
  read(11,*,iostat=ios)
  read(11,*,iostat=ios)
  read(11,*,iostat=ios)
  read(11,FMT='(34x,f6.4)',iostat=ios)propsurf
  read(11,FMT='(26x,f5.2)',iostat=ios)dzeta !previously 25x
  read(11,FMT='(37x,f4.2)',iostat=ios)Bo
  read(11,FMT='(41x,f5.3)',iostat=ios)wid
  bool_rough=1
  if(Bo.eq.0.) then
     opt2='n'
  else
     opt2='y'
  endif
  read(11,*,iostat=ios)
  if(ios/=0) goto 9130

  read(10,FMT='(18x,I1)',iostat=ios)nbl
  read(10,*,iostat=ios)
  if(ios/=0) goto 9115

  allocate(nat(nbl,0:nbcmax),tabcompa(nbl),tabtau(nbl),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='16' 
     goto 9065
  end if
  allocate(tabg(nbl),tabprop(nbl,nbcmax),tabdiam(nbl,nbcmax),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='17' 
     goto 9065
  end if
  allocate(tabs(nbl,nbcmax),tex(nbl),nbpc(nbl),tabthick(nbl),tabwidth(nbl),stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='18' 
     goto 9065
  end if
  if(opt1.eq.'iden') then
     allocate(ap(0:nbord,nbl),stat=err_stat)
     ! check if the allocation is successfull
     if(err_stat/=0) then
        err_code_mem='19' 
        goto 9065
     end if
     allocate(tabomeg(nbl),stat=err_stat)
     ! check if the allocation is successfull
     if(err_stat/=0) then
        err_code_mem='25' 
        goto 9065
     end if
  endif

  nat=0

  loop_layer1: do m=1,nbl

     read(10,*,iostat=ios)
     if(ios/=0) goto 9115
     read(11,*,iostat=ios)
     if(ios/=0) goto 9130

     read(10,FMT='(16x,A1)',iostat=ios)tex(m)
     read(10,FMT='(33x,I2)',iostat=ios)nbpc(m)
     if((nbpc(m)/10).eq.0) then
        nb_digit=1
     else if ((nbpc(m)/10).eq.1) then
        nb_digit=2
     else 
        goto 9117
     endif

     if(nb_digit.eq.1) then
        write(fmt_string,FMT='("(18x,",I1,"(1x,I2))")')nbpc(m)
     else
        write(fmt_string,FMT='("(18x,",I2,"(1x,I2))")')nbpc(m)
     endif

     select case(trim(adjustl(tex(m))))
     case('g')
        read(10,FMT=fmt_string,iostat=ios)(nat(m,1:nbpc(m)))
     case('a')
        read(10,FMT=fmt_string,iostat=ios)(nat(m,1:nbpc(m)))
     case('c')
        read(10,FMT=fmt_string,iostat=ios)(nat(m,0:nbpc(m)-1))
     case default
        goto 9118
     end select
     if(ios/=0) goto 9115

     read(11,FMT='(33x,f10.8)',iostat=ios)tabcompa(m)
     if(m.gt.1) then
        if(opt6.eq.'y') then
           read(11,FMT='(35x,e8.2)',iostat=ios)tabthick(m)
        else
           read(11,FMT='(36x,e8.2)',iostat=ios)tabthick(m)
        endif
        if((opt7.eq.'y').and.(trim(adjustl(tex(m))).eq.'a'))&
             read(11,FMT='(36x,f8.4)',iostat=ios)tabwidth(m)
     else
        tabthick(m)=250.
     endif
     read(11,FMT='(33x,f5.2)',iostat=ios)tabg(m)

     select case(trim(adjustl(tex(m))))
     case('g')
        if(nb_digit.eq.1) then
           write(fmt_string,FMT='("(19x,",I1,"(1x,f14.12))")')nbpc(m)
        else
           write(fmt_string,FMT='("(19x,",I2,"(1x,f14.12))")')nbpc(m)
        endif
        read(11,FMT=fmt_string,iostat=ios)(tabprop(m,1:nbpc(m)))
        if(nb_digit.eq.1) then
           write(fmt_string,FMT='("(23x,",I1,"(1x,e10.5))")')nbpc(m)
        else
           write(fmt_string,FMT='("(23x,",I2,"(1x,e10.5))")')nbpc(m)
        endif
        read(11,FMT=fmt_string,iostat=ios)(tabdiam(m,1:nbpc(m)))
        if(nb_digit.eq.1) then
           write(fmt_string,FMT='("(49x,",I1,"(1x,f6.4))")')nbpc(m)
        else
           write(fmt_string,FMT='("(49x,",I2,"(1x,f6.4))")')nbpc(m)
        endif
        read(11,FMT=fmt_string,iostat=ios)(tabs(m,1:nbpc(m)))
     case('a')
        if(nb_digit.eq.1) then
           write(fmt_string,FMT='("(19x,",I1,"(1x,f14.12))")')nbpc(m)
        else
           write(fmt_string,FMT='("(19x,",I2,"(1x,f14.12))")')nbpc(m)
        endif
        read(11,FMT=fmt_string,iostat=ios)(tabprop(m,1:nbpc(m)))
        if(nb_digit.eq.1) then
           write(fmt_string,FMT='("(23x,",I1,"(1x,e10.5))")')nbpc(m)
        else
           write(fmt_string,FMT='("(23x,",I2,"(1x,e10.5))")')nbpc(m)
        endif
        read(11,FMT=fmt_string,iostat=ios)(tabdiam(m,1:nbpc(m)))
        if(nb_digit.eq.1) then
           write(fmt_string,FMT='("(49x,",I1,"(1x,f6.4))")')nbpc(m)
        else
           write(fmt_string,FMT='("(49x,",I2,"(1x,f6.4))")')nbpc(m)
        endif
        read(11,FMT=fmt_string,iostat=ios)(tabs(m,1:nbpc(m)))
     case('c')
        if(nb_digit.eq.1) then
           write(fmt_string,FMT='("(19x,",I1,"(1x,f14.12))")')nbpc(m)
        else
           write(fmt_string,FMT='("(19x,",I2,"(1x,f14.12))")')nbpc(m)
        endif
        read(11,FMT=fmt_string,iostat=ios)(tabprop(m,1:(nbpc(m)-1)))
        if(nb_digit.eq.1) then
           write(fmt_string,FMT='("(23x,",I1,"(1x,e10.5))")')nbpc(m)
        else
           write(fmt_string,FMT='("(23x,",I2,"(1x,e10.5))")')nbpc(m)
        endif
        read(11,FMT=fmt_string,iostat=ios)(tabdiam(m,1:(nbpc(m)-1)))
        if(nb_digit.eq.1) then
           write(fmt_string,FMT='("(49x,",I1,"(1x,f6.4))")')nbpc(m)
        else
           write(fmt_string,FMT='("(49x,",I2,"(1x,f6.4))")')nbpc(m)
        endif
        read(11,FMT=fmt_string,iostat=ios)(tabs(m,1:(nbpc(m)-1)))
     end select

     if(ios/=0) goto 9130

     read(11,*,iostat=ios)
     if(ios<0) exit
     if(ios>0) goto 9130

     read(10,*,iostat=ios)
     if(ios<0) exit
     if(ios>0) goto 9115

  end do loop_layer1

  !*************************************************************

  !*** Test of the contribution of the current area
  count_nb_area=0
  print *, 'nb wave', nbw
  !print *, 'nb max incid', MAXUMU
  allocate(tab_fact(nbw, MAXUMU), stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='25' 
    goto 9065
  end if
  allocate(tab_sep(nbw, MAXUMU), stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='25' 
    goto 9065
  end if
  allocate(tab_se(nbw, MAXUMU), stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='25' 
    goto 9065
  end if

  allocate(tab_spec(nbw, 32767), stat=err_stat)
  ! check if the allocation is successfull
  if(err_stat/=0) then
     err_code_mem='25' 
    goto 9065
  end if

  dzeta_old=0.
  loop_area: do

     count_nb_k=0
     if(opt1.eq.'diff') then
        write(err_msg,FMT='("  Reading characteristics of area",I5," OK")') nbta
        WRITE(UNIT=6,FMT='(A255)')err_msg
        read(13,*,iostat=ios)
        read(13,FMT='(37x,A1)',iostat=ios)contrib
        if(ios<0) goto 150
        if(ios>0) goto 9135
        read(13,*,iostat=ios)
        if(ios>0) goto 9135
        if(contrib.eq.'n') then
           do 
              read(12,iostat=ios)u,v
              if(ios>0) goto 9140
              if(ios<0) exit                
              if(u==-1) exit
           end do
           write(err_msg,FMT=&
                &'("  Reading spatial distrib. of area",I5," OK")')nbta
           WRITE(UNIT=6,FMT='(A255)')err_msg
           write(err_msg,FMT='("  Area",I4," : not treated")')nbta
           WRITE(UNIT=6,FMT='(A255)')err_msg
           nbta=nbta+1
           goto 150
        endif
     else if(nbta.eq.1) then
        write(err_msg,FMT='("  Reading characteristics of area",I5," OK")')nbta
        WRITE(UNIT=6,FMT='(A255)')err_msg
     endif

     !*************************************************************

     !*** Allocation of the arrays of the Fourier components
     !    and of the single and double scattering corrections

     if(opt1.eq.'diff') then 
        allocate(ap(0:nbord,nbl),stat=err_stat)
        ! check if the allocation is successfull
        if(err_stat/=0) then
           err_code_mem='28' 
           goto 9065
        end if
        allocate(tabomeg(nbl),stat=err_stat)
        ! check if the allocation is successfull
        if(err_stat/=0) then
           err_code_mem='34' 
           goto 9065
        end if
     endif

     !*************************************************************

     !*** Calculation of the Legendre components for the phase
     !    function P

     do m=1,nbl
        g=real(tabg(m),4)
        !   CALL  GETMOM( 3, g, nbord, PMOM(0,nbl-m+1) )
        PMOM(0:nbord,nbl-m+1)=g**expo

     end do


     !*************************************************************

     !*** Initialization of the geographical chain list for 
     !    the current area 

     ! initialization of the geometry counters and arrays
     NUMU0=0
     NUMU=0
     NPHI=0
     NTAU= 0 ! sounding optical depths 
     TUMU0=2.
     TUMU=2.
     TPHI=500.
     TLEVEL = -1.
     mask=.true.
     rang=0

     call init_geo1(geo1ptm%head,geo1ptm%tail)

     allocate(geo1_item,stat=err_stat)
     ! check if the allocation is successfull
     if(err_stat/=0) then
        err_code_mem='37' 
        goto 9065
     end if

     up=-1

     loop_pt_area1: do

        read(12,iostat=ios)u,v
        if(ios>0) goto 9140
        if(ios<0) exit loop_pt_area1 
        if(u==-1) exit loop_pt_area1

        !* change of latitude line in the area geographical description
        if(u/=up) then
           geo2_item=>geo2ptm_tab(u)%head
           vp=1
           up=u
        endif

        !* search of the current geographical point (u,v) in the
        !  latitude line u of the general geographical dynamic 
        !  structure
        do t=1,(v-vp)
           geo2_item=>geo2_item%next
           vp=vp+1
        end do

        muo=geo2_item%muo
        mu=geo2_item%mu
        phi=geo2_item%phi
        alpha=geo2_item%alpha
        depth=geo2_item%depth

        !* only the points both illuminated and observed are treated
        test_valid_pt:if(((muo.gt.0.).and.(muo.le.1.)).and.&
             &((mu.ge.-1.).and.(mu.le.1.)).and.&
             &((alpha.ge.0.).and.(alpha.le.180.)).and.&
             &((phi.ge.(-2.*pi)).and.(phi.le.(2.*pi))).and.&
             &((geo2_item%p.gt.0).and.(geo2_item%q.gt.0))) then

           thetauo=acos(muo)/pi*180.
           thetau=acos(mu)/pi*180.

           !* Calculation of the roughness correction function
           if ((dzeta.gt.0.).and.(dzeta.le.60).and.(mu.gt.0.)) then
              muo_a=muo
              mu_a=mu                 
              call correct(thetauo, thetau, phi, dzeta,& 
                   &muo_a, mu_a, muo_b, mu_b)
              if(((muo_a.eq.1.).or.(mu_a.eq.1.)).or.(alpha.lt.2)) then
                 phi_a=0.
              else
                 phi_a=(-cos(pi*alpha/180.)-muo_a*mu_a)&
                      &/(sqrt(1-muo_a**2)*sqrt(1-mu_a**2))
                 if(phi_a.gt.1.)  phi_a=1.
                 if(phi_a.lt.-1.) phi_a=-1.
                 phi_a=acos(phi_a)
              endif
           else
              muo_a=0.
              mu_a=0.
              muo_b=0.
              mu_b=0.
           endif

           ! modified geometries in case of roughness (Hapke and slab model only)
           geo1_item%muo_a=muo_a
           geo1_item%mu_a=mu_a
           geo1_item%phi_a=phi_a
           geo1_item%muo_b=muo_b
           geo1_item%mu_b=mu_b
           ! untouched geometry
           geo1_item%thetauo=thetauo
           geo1_item%thetau=thetau
           geo1_item%phi=phi
           geo1_item%depth=depth
           geo1_item%alpha=alpha
           ! indexes
           geo1_item%u=u
           geo1_item%v=v
           geo1_item%p=geo2_item%p
           geo1_item%q=geo2_item%q

           if(.not. any(abs(TUMU0-muo).le.1.e-6)) then
              NUMU0=NUMU0+1
              TUMU0(NUMU0)=muo
           endif

           if(.not. any(abs(TUMU-mu).le.1.e-6)) then
              NUMU=NUMU+1
              TUMU(NUMU)=mu
           endif

           if(.not. any(abs(TPHI-phi*180./pi).le.1.e-5)) then
              NPHI=NPHI+1
              TPHI(NPHI)=phi*180./pi
           endif

           if(.not. any(abs(TLEVEL-depth).le.1.e-6)) then
              NTAU=NTAU+1
              TLEVEL(NTAU)=depth
           endif

           call add_geo1(geo1_item,geo1ptm%head,geo1ptm%tail)
           allocate(geo1_item,stat=err_stat)
           ! check if the allocation is successfull
           if(err_stat/=0) then
              err_code_mem='42' 
              goto 9065
           end if

        endif test_valid_pt

     end do loop_pt_area1

     ! sorting the array of emergence angles
     do j=1,NUMU
        tab_rang=minloc(TUMU(1:NUMU),mask(1:NUMU))
        rang(j)=tab_rang(1)
        mask(rang(j))=.false.
     enddo
     TUMUS(1:NUMU)=TUMU(rang(1:NUMU))
     TUMU(1:NUMU)=TUMUS(1:NUMU) 

     if(opt1.eq.'diff') then
        write(err_msg,FMT='("  Reading spatial distrib. of area",I5," OK")')nbta
        WRITE(UNIT=6,FMT='(A255)')err_msg
     else if(nbta.eq.1) then
        write(err_msg,FMT='("  Reading spatial distrib. of area",I5," OK")')nbta
        WRITE(UNIT=6,FMT='(A255)')err_msg
     end if

     !*************************************************************

     ! Reading the file of atmospheric spectral optical depth

     if(opt7.eq.'y') then

        open(15,file=atmtau,form='unformatted',status="old",iostat=ios)
        if(ios/=0) goto 9115

        ! allocation of the optical depth files

        allocate(tauatm(nbw,nbatm+NTAU),stat=err_stat)
        ! check if the allocation is successfull
        if(err_stat/=0) then
           goto 9065
        end if

        do m=1,nbatm+NTAU 
           read(15,iostat=ios)tauatm(:,m)
           if(ios>0) goto 9178
        end do

        close(15)

     endif

     !*************************************************************

     !*************************************************************
     !*************************************************************
     !**                                                         **
     !**  Calculation of the spectral contribution of the        **
     !**  current area to the pixels involved                    **
     !**                                                         **
     !*************************************************************
     !*************************************************************

     !*************************************************************

     !*** Allocation of the array containing the wavelengths
     allocate(tab_w(nbw),stat=err_stat)
     ! check if the allocation is successfull 
     if(err_stat/=0) then
        err_code_mem='51' 
        goto 9065
     end if

     !*** Coming back to the beginnings of the wavelength and 
     !    optical indexes chain lists

     if(.not.associated(wptm%head)) goto 9145     
     w_item=>wptm%head
     do nb=1,nbc
        if(.not.associated(indptm_tab(nb)%head)) goto 9145
        indptm(nb)%head=>indptm_tab(nb)%head
     end do
     nbcw=0

     ! read the current line number
     geo1_item=>geo1ptm%head
     p=geo1_item%p

     perc_mem=0.

     !*************************************************************

     !*** Reading the refractive complex indexes of all the 
     !    compounds for the current wavelength or wave number

     loop_block: do nbc_block=1,nb_block

        sf_nb=((p-1)*nb_block+nbc_block-1)/recpmax+1
        rec_nb=mod((p-1)*nb_block+nbc_block-1,recpmax)+1
        !read(19+sf_nb,REC=rec_nb,iostat=ios)pix_mem
        if(propsurf.ne.1.) then
	  read(19+sf_nb,REC=rec_nb,iostat=ios)pix_mem
          if(ios/=0) goto 9176
        else
          pix_mem=0.
	endif
        nbcwp=1 
        print *,nbwp, nbcw
        !count_nb_k=0

        loop_wavelength: do while(nbcwp.le.nbwp)
   
           wave=w_item%wave 

           ! counting the percentage of calculated wavelengths
           ! for the current area              
           perc_wave=(wave-fw)/(lww-fw)*100.

           perc_comp=perc_wave-mod(perc_wave,10.)

           ! only multiple of 10. are displayed
           if(perc_comp.ne.perc_mem) then
              write(6,fmt='(A23,f5.1,A2)'),'percentage completed : ', perc_comp, ' %'
              perc_mem=perc_comp
	      !print *, 'nb sep sampled : ', count_nb_k
           endif


        
           do nb=1,nbc
              ind_item=>indptm(nb)%head
              tabindreal(nb)=ind_item%indreal
              tabindimag(nb)=ind_item%indimag
           end do

           !!!!!!!!!!!!!!!!!!!!
           !Calculation of dn and dk used to determine if we need to recalculate sep 
           !at the wavelenght nbcwp (modif F.A)
           if(nbcwp.eq.1) then
              last_tabindreal=tabindreal(nat(2,0))  !initialization
              last_tabindimag=tabindimag(nat(2,0))
           endif
           tabn(nbcwp)=tabindreal(nat(2,0))
           tabk(nbcwp)=tabindimag(nat(2,0))
           dn= abs( (tabn(nbcwp)-last_tabindreal)/last_tabindreal )
           dk=abs( (log(tabk(nbcwp))-log(last_tabindimag))/log(last_tabindimag) )
           last_tabindreal=tabn(nbcwp) !store the last value of tabindreal
           last_tabindimag=tabk(nbcwp) 

           !!!!!!!!!!!!!!!!!!!!!!


           nbcw=nbcw+1
           if(nbcw.eq.nbw) lw=wave

           ! initializing the user optical depths 
           UTAU=0.
           if(opt7.eq.'y') then
              altitude=0. ! the altitude level
              sumtau=0. ! the total optical depth (particles + aerosols)
              do j=1,NTAU
                 ind_lyr(j)=nbl-nbatm
                 tab_alt(j)=0.
              end do
           endif

           !*************************************************************

           !*** Calculation of the radiative transfer quantities 
           !    independent to the local conditions of illumination
           !    and observation for the current area

           loop_layer2: do m=1,nbl

              !* calculation of the single scattering albedo and the 
              !  optical depth for the current layer

              if(tex(m).eq.'g') then
                 !  granular texture

                 allocate(indr(nbpc(m)),indi(nbpc(m)),prop(nbpc(m)),&
                      &stat=err_stat)
                 ! check if the allocation is successfull
                 if(err_stat/=0) then
                    err_code_mem='43' 
                    goto 9065
                 end if
                 allocate(diam(nbpc(m)),s(nbpc(m)),stat=err_stat)
                 ! check if the allocation is successfull
                 if(err_stat/=0) then
                    err_code_mem='44' 
                    goto 9065
                 end if

                 indr=tabindreal(nat(m,1:nbpc(m)))
                 indi=tabindimag(nat(m,1:nbpc(m)))
                 prop=tabprop(m,1:nbpc(m))
                 diam=tabdiam(m,1:nbpc(m))
                 s=tabs(m,1:nbpc(m))

                 if((dim1.eq.'cm').and.(m.eq.1)) wave=1./wave*1.d4
                 if((dim1.eq.'nm').and.(m.eq.1)) wave=wave/1.d3

                 call proptexg(nbpc(m),indr,indi,prop,diam,s,&
                      &tabcompa(m),wave,tabthick(m),singalbg,tau)
                 tabomeg(nbl-m+1)=singalbg

                 if(opt6.eq.'y') then
                    if(m==1) then
                       tabtau(nbl)=250.
                    else
                       tabtau(nbl-m+1)=tabthick(m)
                    endif
                 else
                    if(m==1) then
                       tabtau(nbl)=250.
                    else
                       tabtau(nbl-m+1)=tau
                    endif
                 endif

                 deallocate(indr,indi,prop,diam,s)

              else if (tex(m).eq.'c') then
                 !  compact texture

                 allocate(indr(0:(nbpc(m)-1)),indi(0:(nbpc(m)-1)),&
                      &stat=err_stat)
                 ! check if the allocation is successfull
                 if(err_stat/=0) then
                    err_code_mem='45' 
                    goto 9065
                 end if
                 allocate(prop(nbpc(m)-1),diam(nbpc(m)-1),stat=err_stat)
                 ! check if the allocation is successfull
                 if(err_stat/=0) then
                    err_code_mem='46' 
                    goto 9065
                 end if
                 allocate(s(nbpc(m)-1),stat=err_stat)
                 ! check if the allocation is successfull
                 if(err_stat/=0) then
                    err_code_mem='47' 
                    goto 9065
                 end if

                 indr=tabindreal(nat(m,0:nbpc(m)-1))
                 indi=tabindimag(nat(m,0:nbpc(m)-1))
                 prop=tabprop(m,1:nbpc(m)-1)
                 diam=tabdiam(m,1:nbpc(m)-1)
                 s=tabs(m,1:nbpc(m)-1)

                 if((dim1.eq.'cm').and.(m.eq.1)) wave=1./wave*1.d4
                 if((dim1.eq.'nm').and.(m.eq.1)) wave=wave/1.d3

                 call proptexc(int2(nbpc(m)-1),indr,indi,prop,diam,s,&
                      &tabcompa(m),wave,tabthick(m),singalbc,tau, R0, T0)

                 tabomeg(nbl-m+1)=singalbc
                 if(m==1) then
                    tabtau(nbl)=250.
                 else
                    tabtau(nbl-m+1)=tau
                 endif

                 deallocate(indr,indi,prop,diam,s)

              else if (tex(m).eq.'a') then
                 ! atmospheric layer

                 if(tabthick(m).gt.0.) then

                    allocate(indr(nbpc(m)),indi(nbpc(m)),prop(nbpc(m)),&
                         &stat=err_stat)
                    ! check if the allocation is successfull
                    if(err_stat/=0) then
                       err_code_mem='43' 
                       goto 9065
                    end if
                    allocate(diam(nbpc(m)),s(nbpc(m)),stat=err_stat)
                    ! check if the allocation is successfull
                    if(err_stat/=0) then
                       err_code_mem='44' 
                       goto 9065
                    end if

                    indr=tabindreal(nat(m,1:nbpc(m)))
                    indi=tabindimag(nat(m,1:nbpc(m)))
                    prop=tabprop(m,1:nbpc(m))
                    diam=tabdiam(m,1:nbpc(m))
                    s=tabs(m,1:nbpc(m))

                    if((dim1.eq.'cm').and.(m.eq.1)) wave=1./wave*1.d4
                    if((dim1.eq.'nm').and.(m.eq.1)) wave=wave/1.d3

                    call proptexg(nbpc(m),indr,indi,prop,diam,s,&
                         &tabcompa(m),wave,tabthick(m),singalbg,tau)

                    if(opt6.eq.'y') then

                       tabtau(nbl-m+1)=tabthick(m)+tauatm(nbcw,m-nbatm+1)
                       tau=tabthick(m)

                    else

                       tabtau(nbl-m+1)=tau+tauatm(nbcw,m-nbatm+1)

                    endif

                    altitude=altitude+tabwidth(m)
                    ! summing the optical depth of the particles to calculate the user optical depths
                    do j=1,NTAU
                       if(TLEVEL(j).GE.altitude) then
                          UTAU(j)=UTAU(j)+tau
                          ind_lyr(j)=m
                          tab_alt(j)=altitude
                       endif
                    end do

                    betap=tau/tabwidth(m)
                    betag=tauatm(nbcw,m-nbatm+1)/tabwidth(m)
                    fact_D=betag/betap
                    singalbg=singalbg/(1+fact_D)
                    deallocate(indr,indi,prop,diam,s)                       
                 else
                    tabtau(nbl-m+1)=tauatm(nbcw,m-nbatm+1)
                    singalbg=0.
                 endif
                 tabomeg(nbl-m+1)=singalbg

              endif

           enddo loop_layer2

           ! calculation of the total optical depth of the atmosphere
           ! (particles+gas) 

           if(opt7.eq.'y') then 
              do m=1,nbatm              
                 sumtau=sumtau+tabtau(m)                 
              end do

              do j=1,NTAU
                 ! summing the optical depth of the gas to calculate the user optical depths
                 UTAU(j)=UTAU(j)+tauatm(nbcw,nbatm+j)
                 ! and the remaining optical depth linked to the particles
                 if(tabthick(ind_lyr(j)+1).gt.0.) UTAU(j)=UTAU(j)+tabthick(ind_lyr(j)+1)*&
                      &(TLEVEL(j)-tab_alt(j))/tabwidth(ind_lyr(j)+1)
                 ! inversing the 0 optical depth level (ground -> top of the atmosphere)
                 UTAU(j)=sumtau-UTAU(j)
                 if (UTAU(j).lt.0.) then 
                    UTAU(j)=0. 
                 else if (UTAU(j).gt.sumtau) then 
                    UTAU(j)=sumtau
                 endif

              end do

           endif


           !*************************************************************

           !*** Integration of the area spectral contribution points 
           !    by points

           if(.not.associated(geo1ptm%head)) cycle loop_area
           geo1_item=>geo1ptm%head

           loop_inc: do i=1,NUMU0

              FBEAM= pi/TUMU0(i)

              ! *** DISORT algorithm
		!MODIF F.Andrieu 08/2013 condition for disort running nbl>2 --> nbl>3
              if((nbl>3).or.(opt5=='n')) then

                 CALL  DISORT( nbl, real(tabtau,4), real(tabomeg,4), nbord, PMOM, TEMPER,&
                      & WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,&
                      & USRANG, NUMU, TUMU, NPHI, TPHI, IBCND, FBEAM,&
                      & TUMU0(i), PHI0, FISOT, LAMBER, ALBEDO, BTEMP,&
                      & TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,&
                      & HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,&
                      & MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,&
                      & ALBMED, TRNMED )

              endif

              ! *** slab over granular substrat model. The slab can be covered by a thin granular layer.

              if(((nbl==2).or.(nbl==3)).and.(opt5=='y').and.(tex(2).eq.'c')) then
 

                 indrs=real(tabindreal(nat(2,0)),4) 
                 indis=real(tabindimag(nat(2,0)),4) 
                 theta0=acos(TUMU0(i))*180./pi   
		

           !Modif F.ANDRIEU 07/2014 Changing the way low reolution calculation are done, 
           !avoiding missing spectral bans (slighty longer)

                 if(bool_rough.eq.1) then 
		   if(  (nbcwp.eq.1).or.(opt4.eq.'d').or.(dn>0.01).or.(dk>0.1)  ) then 
                    !CO2 :  (dn>0.01).or.(dk>0.1) ; Water : (dn>0.000075).or.(dk>0.00075)
                    count_nb_k=count_nb_k+1
                    fact=est_fact(indrs,theta0,dzeta)
                    tab_fact(nbcwp+(nbc_block-1)*nbwp, i)=fact
                    sep=coefreflgp(indrs, indis, theta0, dzeta)
		    if(sep.gt.1.) then
                      sep=1.
                      if(nbcwp.eq.1) print *, 'Error in Sep estimation for incidence (°)', theta0
                    endif
                    tab_sep(nbcwp+(nbc_block-1)*nbwp, i)=sep
                    call coefreflcp(tabindreal(nat(2,0)), tabindimag(nat(2,0)), se, si)
                    tab_se(nbcwp+(nbc_block-1)*nbwp, i)=se 
		  else
		    fact=tab_fact(nbcwp+(nbc_block-1)*nbwp-1, i)
                    sep=tab_sep(nbcwp+(nbc_block-1)*nbwp-1, i)
                    se=tab_se(nbcwp+(nbc_block-1)*nbwp-1, i) 
		    tab_fact(nbcwp+(nbc_block-1)*nbwp, i)=fact
		    tab_sep(nbcwp+(nbc_block-1)*nbwp, i)=sep
		    tab_se(nbcwp+(nbc_block-1)*nbwp, i)=se
		   endif
                 else
		    !if a new estimation is not necessary, we keep the last estimation
		    !N.B. : index i belong to the solar incidence loop

                    fact=tab_fact(nbcwp+(nbc_block-1)*nbwp, i)
                    sep=tab_sep(nbcwp+(nbc_block-1)*nbwp, i)
                    se=tab_se(nbcwp+(nbc_block-1)*nbwp, i) 
                 endif

                 call slab_refl_trans(real(wave,4), indrs, indis, &
                      & real(tabomeg(nbl-1),4), real(tabtau(nbl-1),4), &
                      &  real(tabthick(2),4), fact, sep, R0p, T0p)

                 R_diff=fct_reflectance_SC_diff(real(R0,4),real(T0,4), &
                      & R0p,T0p,real(tabomeg(nbl),4)) !Total reflexion coefficient
                 if(NTAU.eq.2) then
                    if(nbl==2) then
                       rc=0.
                       tc=1.
                    else
                       gamma=sqrt(1.-real(tabomeg(1),4))
                       ru=(1-gamma)/(1+gamma)
                       rc=ru*(1-exp(-4.*gamma*real(tabtau(1),4)))/ &
                            (1-ru**2*exp(-4.*gamma*real(tabtau(1),4)))
                       tc=(1-ru**2)*exp(-2.*gamma*real(tabtau(1),4))/ &
                            (1-ru**2*exp(-4.*gamma*real(tabtau(1),4)))
                    endif
                    gamma=sqrt(1.-real(tabomeg(nbl),4))
                    rs=(1-gamma)/(1+gamma)
                    T_diff=TUMU0(i)*tc*T0p/(1-real(R0,4)*(rc+rs)+rs*rc*(real(R0,4)**2-real(T0,4)**2))

                 endif

              endif

              loop_depth: do r=1,NTAU 

                 loop_az: do t=1,NPHI

                    loop_emerg: do j=1,NUMU

                       ! *** DISORT algorithm
                       if((nbl>1).or.(opt5=='n')) then

                          if (t==1) then
                             refl=real(UU(rang(j),r,t),8)
                          else
                             refl=real(UU(rang(NUMU-j+1),r,t),8)
                          endif

                       endif

                       muo_a=geo1_item%muo_a
                       mu_a=geo1_item%mu_a
                       muo_b=geo1_item%muo_b
                       mu_b=geo1_item%mu_b
                       phi_a=geo1_item%phi_a
                       thetauo=geo1_item%thetauo
                       thetau=geo1_item%thetau           
                       phi=geo1_item%phi
                       depth=geo1_item%depth
                       alpha=geo1_item%alpha
                       u=geo1_item%u
                       p=geo1_item%p
                       q=geo1_item%q

                       ! *** Hapke Model
                       if((nbl==1).and.(opt5=='y').and.(mu.gt.0.).and.(depth.eq.0.)) then

                          if (spherical.eq.1) then
                             omeg=tabomeg(1)
                             refl=( 1.-sqrt(1.-omeg) )/( 1.+sqrt(1.-omeg) )
                          else

                          g=tabg(nbl)

                          ! Correction of the roughness effect if necessary
                          if((dzeta.gt.0.).and.(dzeta.le.60)) then
                             theta0=acos(muo_a)*180./pi
                             theta=acos(mu_a)*180./pi    
                             refl=real(fct_reflectance(real(tabomeg(1),4), g, 0., theta0, &
                                  &theta, (180.-alpha), err_flag),8)                            
                             refl=refl*real(&
                                  &fct_S(thetauo, thetau, muo_a, mu_a, muo_b, mu_b, phi, dzeta),8)
                          else
                             refl=real(fct_reflectance(real(tabomeg(1),4), g, 0., thetauo, &
                                  &thetau, (180.-alpha), err_flag),8)                                
                          endif


                          !* adding the correction of the single scattering contribution
                          !  for the reflected part by the last layer (opposition effect : eq 2.8-1)
                          if(opt2.eq.'y') then
                             omeg=tabomeg(1)
                             g=tabg(1)
                             refl=refl+(Bo*(1.-tan(pi-rad(alpha))/(2.*wid)&
                                  &*(3.-exp(-wid/tan(pi-rad(alpha))))&
                                  &*(1.-exp(-wid/tan(pi-rad(alpha))))))&
                                  &*fctP(dble(g),dble(rad(alpha)))&
                                  &*omeg/(mu+muo)/4.

                          endif
                          
                          endif
                       endif

                       ! *** slab over granular substrat model
  

                       if((nbl==2).and.(opt5=='y').and.(mu.gt.0.) &
                            .and.(tex(2).eq.'c').and.(depth.eq.0.)) then

                          if(r==1) then
			     if((dzeta.gt.0.).and.(dzeta.le.60)) then
                             	theta0=acos(muo_a)*180./pi
                             	theta=acos(mu_a)*180./pi
                                if(bool_rough.eq.1) then 
		                  if(  (nbcwp.eq.1).or.(opt4.eq.'d').or.(dn>0.01).or.(dk>0.1) ) then 
                             	    R_spec=fct_reflectance_SC_spec(indrs, indis, dzeta, &
                                      & thetauo, thetau, phi, (180.-alpha), err_flag) 
                                    tab_spec(nbcwp+(nbc_block-1)*nbwp, j+NUMU*(t-1)+NUMU*NPHI*(i-1) )=R_spec
                                  else
                                    R_spec=tab_spec(nbcwp+(nbc_block-1)*nbwp-1, j+NUMU*(t-1)+NUMU*NPHI*(i-1) ) 
                                    tab_spec(nbcwp+(nbc_block-1)*nbwp, j+NUMU*(t-1)+NUMU*NPHI*(i-1) )=R_spec  
                                  endif
                                else
                                  R_spec=tab_spec(nbcwp+(nbc_block-1)*nbwp, j+NUMU*(t-1)+NUMU*NPHI*(i-1) )
                                endif
		                refl=real(R_spec+R_diff,8)
			     else
			        refl=real(R_diff,8)
			     endif
                          else
                             refl=real(T_diff,8)
                          endif


                       endif

                       ! *** slab over granular substrat model. The slab itself is covered by a thin granular layer.

                       if((nbl==3).and.(opt5=='y').and.(mu.gt.0.) &
                            .and.(tex(2).eq.'c').and.(depth.eq.0.)) then

                          if(r==1) then
                             g=tabg(nbl)
			     if((dzeta.gt.0.).and.(dzeta.le.60)) then
                             	theta0=acos(muo_a)*180./pi
                             	theta=acos(mu_a)*180./pi
                             	ru=fct_reflectance(real(tabomeg(1),4), g, 0., theta0, &
                                  &theta, (180.-alpha), err_flag) 

                             	gamma=sqrt(1.-real(tabomeg(1),4))

                             	rl=R_diff+(sep+real(se,4))/2.
                             	refl=real(ru*((1.+1./ru*(rl-ru)/(1.-ru*rl)* &
                                  &   exp(-4.*gamma*real(tabtau(1),4)))/ &
                                  & (1.+ru*(rl-ru)/(1.-ru*rl)*exp(-4.*gamma*real(tabtau(1),4)))),8)
				refl=refl*real(&
                                  &fct_S(thetauo, thetau, muo_a, mu_a, muo_b, mu_b, phi, dzeta),8)
			     else
				ru=fct_reflectance(real(tabomeg(1),4), g, 0., thetauo, &
                                  &thetau, (180.-alpha), err_flag) 

                             	gamma=sqrt(1.-real(tabomeg(1),4))

                             	rl=R_diff+(sep+real(se,4))/2.
                             	refl=real(ru*((1.+1./ru*(rl-ru)/(1.-ru*rl)* &
                                  &   exp(-4.*gamma*real(tabtau(1),4)))/ &
                                  & (1.+ru*(rl-ru)/(1.-ru*rl)*exp(-4.*gamma*real(tabtau(1),4)))),8)
			     endif
                          else
                             refl=real(T_diff,8)
                          endif

			 !* adding the correction of the single scattering contribution
                          !  for the reflected part by the last layer
                          if(opt2.eq.'y') then
                             omeg=tabomeg(1)
                             g=tabg(1)
                             refl=refl+(Bo*(1.-tan(pi-rad(alpha))/(2.*wid)&
                                  &*(3.-exp(-wid/tan(pi-rad(alpha))))&
                                  &*(1.-exp(-wid/tan(pi-rad(alpha))))))&
                                  &*fctP(dble(g),dble(rad(alpha)))&
                                  &*omeg/(mu+muo)/4.

                          endif

                       endif

                       !*************************************************************

                       !*** Calculation of the spectral contribution of the current 
                       !    geographical point, for the current wavelength or 
                       !    wavenumber to the corresponding pixels

                       !* calculation of the surface attached to the geographical 
                       !  point

                       !* integration of the contribution of the current surface
                       pix(q,p)=pix(q,p)+refl

                       !* integration of the contribution of a reference lambertian 
                       !  corresponding surface
                       if(nbcw.eq.1)& 
                            &pixa(q,p)=1.

                       if(.not.associated(geo1_item%next)) exit loop_inc

                       geo1_item=>geo1_item%next


                    end do loop_emerg
                 end do loop_az
              end do loop_depth
           end do loop_inc

           !*************************************************************

           !*** Adding the spectral contribution of the current area 
           !    to the pixels involved

           do q=1,qmax
              if(.not.(pix(q,p).eq.0.)) then
                 !* adding the spectral contribution to the current spatial-spectral plan
                 pix_mem(nbcwp,q)=pix_mem(nbcwp,q)+real(pix(q,p)*propsurf)
                 pix(q,p)=0.
              endif
           end do

           !Modif F.ANDRIEU 07/2014 Changing the way low reolution calculation are done, 
           !avoiding missing spectral bans (slighty longer)

!!!!!!!!!! Old version based on the derivative
!           if((nbcwp.gt.1).and.(opt4.ne.'d')) then
!              delta_R_moy=0.
!              do q=1,qmax
!                 delta_R_moy=delta_R_moy+(pix_mem(nbcwp,q)-pix_mem(nbcwp-1,q))
!              end do
!              delta_R_moy=abs(delta_R_moy/qmax)
!              if(delta_R_moy.le.delta_R_max) then
!                 nb_pt_skip=int((1-nb_pt_skip_max)/delta_R_max* &
!                      & delta_R_moy+nb_pt_skip_max)
!                 if(nb_pt_skip.lt.1) nb_pt_skip=1
!              else
!                 nb_pt_skip=1
!              endif
!           else
!              nb_pt_skip=1
!           endif
!
!!!!!!!!!!!New version based on precalculated thresholds (must be conducted for every new material)
	   nb_pt_skip=1
           do n_pt_skip=0,nb_pt_skip-1
              if(.not.associated(w_item%next)) exit loop_wavelength
              w_item=>w_item%next

              do nb=1,nbc
                 ind_item=>indptm(nb)%head
                 if(.not.associated(ind_item%next)) exit loop_wavelength
                 indptm(nb)%head=>ind_item%next
              end do
              do q=1,qmax
                 pix_mem(nbcwp+n_pt_skip,q)=pix_mem(nbcwp,q)
              end do
           end do
!!!!!!!!!!!
           nbcwp=nbcwp+nb_pt_skip

        enddo loop_wavelength
        print *,'nb sep sampled = ', count_nb_k, sep
        sf_nb=((p-1)*nb_block+nbc_block-1)/recpmax+1
        rec_nb=mod((p-1)*nb_block+nbc_block-1,recpmax)+1
        write(19+sf_nb,REC=rec_nb,iostat=ios)pix_mem
        if(ios/=0) goto 9176

     end do loop_block

     !*************************************************************

     !*** Deallocation of the limit arrays
     !    and of the geographical chain of the area treated

     if(opt1.eq.'diff') then 
        deallocate(ap)
        deallocate(tabomeg)
     endif

     do
        call delete_geo1(geo1ptm%head,geo1ptm%tail,geo1_item)
        nullify(geo1_item)
        if(.not.associated(geo1ptm%head)) exit
     enddo

     write(err_msg,FMT='(" Area",I4," : treated")') nbta
     WRITE(UNIT=6,FMT='(A255)')err_msg
     !print *, 'nb sep sampled = ', count_nb_k
     nbta=nbta+1

     !*************************************************************

     !*** Deallocation of the characteristic arrays

150  if(opt1.eq.'diff') then
        deallocate(nat,tabcompa,tabtau,tabg,tabprop)
        deallocate(tabdiam,tabs,tex,nbpc,tabthick, tabwidth)
     endif

     if(opt7.eq.'y') deallocate(tauatm)

     !*************************************************************

     !*************************************************************
     !*************************************************************
     !**                                                         **
     !**  Preparation of the data related to the following area  **
     !**                                                         **
     !*************************************************************
     !*************************************************************

     !*** Initialization of the layer characteristic arrays
     !    for the following area

     if(opt1.eq.'diff') then
        read(10,*,iostat=ios)
        read(10,*,iostat=ios)
        read(10,*,iostat=ios)
        if(ios>0) goto 9115
        if(ios<0) exit loop_area
     endif

     read(11,*,iostat=ios)
     read(11,*,iostat=ios)
     read(11,*,iostat=ios)
     read(11,FMT='(34x,f6.4)',iostat=ios)propsurf
     dzeta_old=dzeta
     read(11,FMT='(26x,f5.2)',iostat=ios)dzeta !Warning : previously 25x, contradictory with create_sfiles
     if (dzeta_old.ne.dzeta) then
       bool_rough=1
     else
       bool_rough=0
     endif
     read(11,FMT='(37x,f4.2)',iostat=ios)Bo
     read(11,FMT='(41x,f5.3)',iostat=ios)wid
     if(Bo.eq.0.) then
        opt2='n'
     else
        opt2='y'
     endif
     read(11,*,iostat=ios)
     if(ios>0) goto 9130
     if(ios<0) exit loop_area

     if(opt1.eq.'diff') then

        read(10,FMT='(18x,I1)',iostat=ios)nbl
        read(10,*,iostat=ios)
        if(ios>0) goto 9115
        if(ios<0) exit loop_area

        allocate(nat(nbl,0:nbcmax),tabcompa(nbl),tabtau(nbl),&
             &stat=err_stat)
        ! check if the allocation is successfull
        if(err_stat/=0) then
           err_code_mem='48' 
           goto 9065
        end if
        allocate(tabg(nbl),tabprop(nbl,nbcmax),tabdiam(nbl,nbcmax),&
             &stat=err_stat)
        ! check if the allocation is successfull
        if(err_stat/=0) then
           err_code_mem='49' 
           goto 9065
        end if
        allocate(tabs(nbl,nbcmax),tex(nbl),nbpc(nbl),tabthick(nbl),&
             &tabwidth(nbl),stat=err_stat)
        ! check if the allocation is successfull
        if(err_stat/=0) then
           err_code_mem='50' 
           goto 9065
        end if

        nat=0

     endif

     loop_layer3: do m=1,nbl

        if(opt1.eq.'diff') then
           read(10,*,iostat=ios)
           if(ios/=0) goto 9115
           read(10,FMT='(16x,A1)',iostat=ios)tex(m)
           read(10,FMT='(33x,I2)',iostat=ios)nbpc(m)
           if((nbpc(m)/10).eq.0) then
              nb_digit=1
           else if ((nbpc(m)/10).eq.1) then
              nb_digit=2
           else 
              goto 9117
           endif

           if(nb_digit.eq.1) then
              write(fmt_string,FMT='("(18x,",I1,"(1x,I2))")')nbpc(m)
           else
              write(fmt_string,FMT='("(18x,",I2,"(1x,I2))")')nbpc(m)
           endif

           select case(trim(adjustl(tex(m))))
           case('g')
              read(10,FMT=fmt_string,iostat=ios)(nat(m,1:nbpc(m)))
           case('a')
              read(10,FMT=fmt_string,iostat=ios)(nat(m,1:nbpc(m)))
           case('c')
              read(10,FMT=fmt_string,iostat=ios)(nat(m,0:nbpc(m)-1))
           case default
              goto 9118
           end select
           if(ios/=0) goto 9115
        endif
        read(11,*,iostat=ios)
        if(ios/=0) goto 9130

        read(11,FMT='(33x,f10.8)',iostat=ios)tabcompa(m)
        if(m.gt.1) then
           if(opt6.eq.'y') then
              read(11,FMT='(35x,e8.2)',iostat=ios)tabthick(m)
           else
              read(11,FMT='(36x,e8.2)',iostat=ios)tabthick(m)
           endif
           if((opt7.eq.'y').and.(trim(adjustl(tex(m))).eq.'a'))&
                read(11,FMT='(36x,f8.4)',iostat=ios)tabwidth(m)
        else
           tabthick(m)=250.
        endif
        read(11,FMT='(33x,f5.2)',iostat=ios)tabg(m)

        select case(trim(adjustl(tex(m))))
        case('g')
           if(nb_digit.eq.1) then
              write(fmt_string,FMT='("(19x,",I1,"(1x,f14.12))")')nbpc(m)
           else
              write(fmt_string,FMT='("(19x,",I2,"(1x,f14.12))")')nbpc(m)
           endif
           read(11,FMT=fmt_string,iostat=ios)&
                &(tabprop(m,1:nbpc(m)))
           if(nb_digit.eq.1) then
              write(fmt_string,FMT='("(23x,",I1,"(1x,e10.5))")')nbpc(m)
           else
              write(fmt_string,FMT='("(23x,",I2,"(1x,e10.5))")')nbpc(m)
           endif
           read(11,FMT=fmt_string,iostat=ios)&
                &(tabdiam(m,1:nbpc(m)))
           if(nb_digit.eq.1) then
              write(fmt_string,FMT='("(49x,",I1,"(1x,f6.4))")')nbpc(m)
           else
              write(fmt_string,FMT='("(49x,",I2,"(1x,f6.4))")')nbpc(m)
           endif
           read(11,FMT=fmt_string,iostat=ios)&
                &(tabs(m,1:nbpc(m)))
        case('c')
           if(nb_digit.eq.1) then
              write(fmt_string,FMT='("(19x,",I1,"(1x,f14.12))")')nbpc(m)
           else
              write(fmt_string,FMT='("(19x,",I2,"(1x,f14.12))")')nbpc(m)
           endif
           read(11,FMT=fmt_string,iostat=ios)&
                &(tabprop(m,1:(nbpc(m)-1)))
           if(nb_digit.eq.1) then
              write(fmt_string,FMT='("(23x,",I1,"(1x,e10.5))")')nbpc(m)
           else
              write(fmt_string,FMT='("(23x,",I2,"(1x,e10.5))")')nbpc(m)
           endif
           read(11,FMT=fmt_string,iostat=ios)&
                &(tabdiam(m,1:(nbpc(m)-1)))
           if(nb_digit.eq.1) then
              write(fmt_string,FMT='("(49x,",I1,"(1x,f6.4))")')nbpc(m)
           else
              write(fmt_string,FMT='("(49x,",I2,"(1x,f6.4))")')nbpc(m)
           endif
           read(11,FMT=fmt_string,iostat=ios)&
                &(tabs(m,1:(nbpc(m)-1)))
        case('a')
           if(nb_digit.eq.1) then
              write(fmt_string,FMT='("(19x,",I1,"(1x,f14.12))")')nbpc(m)
           else
              write(fmt_string,FMT='("(19x,",I2,"(1x,f14.12))")')nbpc(m)
           endif
           read(11,FMT=fmt_string,iostat=ios)&
                &(tabprop(m,1:nbpc(m)))
           if(nb_digit.eq.1) then
              write(fmt_string,FMT='("(23x,",I1,"(1x,e10.5))")')nbpc(m)
           else
              write(fmt_string,FMT='("(23x,",I2,"(1x,e10.5))")')nbpc(m)
           endif
           read(11,FMT=fmt_string,iostat=ios)&
                &(tabdiam(m,1:nbpc(m)))
           if(nb_digit.eq.1) then
              write(fmt_string,FMT='("(49x,",I1,"(1x,f6.4))")')nbpc(m)
           else
              write(fmt_string,FMT='("(49x,",I2,"(1x,f6.4))")')nbpc(m)
           endif
           read(11,FMT=fmt_string,iostat=ios)&
                &(tabs(m,1:nbpc(m)))
        case default
           goto 9118
        end select

        if(ios/=0) goto 9130

        if(opt1.eq.'diff') then
           read(10,*,iostat=ios)
           if(ios>0) goto 9115
           if(ios<0) cycle loop_area
        endif

        read(11,*,iostat=ios)
        if(ios>0) goto 9130
        if(ios<0) cycle loop_area

     end do loop_layer3
  
     count_nb_area=count_nb_area+1
 
  end do loop_area

  !*************************************************************

  !*************************************************************
  !*************************************************************
  !**                                                         **
  !**  Destruction of the unused data before the last         **
  !**  treatments                                             **
  !**                                                         **
  !*************************************************************
  !*************************************************************

  !*** Deallocation of the data arrays and the data chain lists

  if(opt1.eq.'iden') then
     deallocate(nat,tabcompa,tabtau,tabg,tabprop)
     deallocate(tabdiam,tabs,tex,nbpc,tabthick, tabwidth)
     deallocate(ap)
     deallocate(tabomeg)
  endif

  do u=1,umax
     do
        call delete_geo2(geo2ptm_tab(u)%head,geo2ptm_tab(u)%tail,geo2_item)
        nullify(geo2_item)
        if(.not.associated(geo2ptm_tab(u)%head)) exit
     enddo
  end do

  deallocate(geo2ptm_tab)    

  deallocate(tabindreal,tabindimag,tabdindreal,tabdindimag) !,last_tabindreal,last_tabindimag
  deallocate(pix)

  do nb=1,nbc
     do
        call delete_ind(indptm_tab(nb)%head,indptm_tab(nb)%tail,ind_item)
        nullify(ind_item)
        if(.not.associated(indptm_tab(nb)%head)) exit
     enddo
  end do

  deallocate(indptm_tab,indptm)

  deallocate(tab_fact)
  deallocate(tab_sep)
  deallocate(tab_se)

  !*************************************************************

  !*** Closing some data files

  close(10)
  close(11)
  close(12)
  close(13)

  !*************************************************************
  !*** 

  write(err_msg,FMT='("**** Ending the calculations ****")')
  WRITE(UNIT=6,FMT='(A255)')err_msg

  write(err_msg,FMT='("**** Pushing the result into the raw binary file ****")')
  WRITE(UNIT=6,FMT='(A255)')err_msg

     ! destroying the wavelength pile
  do nbcw=1,nbw
     call delete_w(wptm%head,wptm%tail,w_item)
     tab_w(nbcw)=w_item%wave
     nullify(w_item)
     if(.not.associated(wptm%head)) exit
  enddo

  write(UNIT=19,iostat=ios)tab_w
  if(ios/=0) goto 9176

  write(UNIT=19,iostat=ios)pixa
  if(ios/=0) goto 9176

  if(opt3.eq.'raw') then 
       
     close(14)

     ! create the raw binary output file
     ! opening
     open(unit=30, file=raw_output, form='unformatted', status='new', iostat=ios)
     if(ios/=0) goto 9179
     ! writing header
     write(unit=30)pmax,qmax,nbw
     write(unit=30)tab_w

     ! writing raw spectra 
     if(nb_block.gt.1) goto 9180
     nbc_block=1

     do p=1,pmax
        sf_nb=((p-1)*nb_block+nbc_block-1)/recpmax+1
        rec_nb=mod((p-1)*nb_block+nbc_block-1,recpmax)+1
        read(19+sf_nb,REC=rec_nb,iostat=ios)pix_mem
        if(ios/=0) goto 9176  
        do q=1,qmax
           if(pixa(q,p)==0.) then 
              write(unit=30)-1.
           else            
              select case(dim1)
              case('mi')
                 select case(dim2)
                 case('mi')
                    write(unit=30)pix_mem(:,q)/real(pixa(q,p))
                 case('cm')
                    write(unit=30)(pix_mem(nbw-nbcw+1,q)/real(pixa(q,p)),nbcw=1,nbw) 
                 case('nm')
                    write(unit=30)pix_mem(:,q)/real(pixa(q,p))
                 end select
              case('cm')
                 select case(dim2)
                 case('mi')
                    write(unit=30)(pix_mem(nbw-nbcw+1,q)/real(pixa(q,p)),nbcw=1,nbw)
                 case('cm')
                    write(unit=30)pix_mem(:,q)/real(pixa(q,p))
                 case('nm')
                    write(unit=30)(pix_mem(nbw-nbcw+1,q)/real(pixa(q,p)),nbcw=1,nbw)
                 end select
              case('nm')
                 select case(dim2)
                 case('mi')
                    write(unit=30)pix_mem(:,q)/real(pixa(q,p))
                 case('cm')
                    write(unit=30)(pix_mem(nbw-nbcw+1,q)/real(pixa(q,p)),nbcw=1,nbw)
                 case('nm')
                    write(unit=30)pix_mem(:,q)/real(pixa(q,p))
                 end select
              end select
           endif
        enddo
     enddo

     ! closing
     close(unit=30)

  endif

  write(err_msg,FMT='("******<< Ending the simulation >>******")')
  WRITE(UNIT=6,FMT='(A255)')err_msg

  !*************************************************************

   !*** Deallocation of the data arrays
  deallocate(pixa,pix_mem,pix_bck)
  deallocate(tab_w)

   !*** Closing the scratch files
  do k=0,nb_sfiles
     close(19+k)
  end do

stop

!*************************************************************
!*************************************************************
!**                                                         **
!**  Error messages                                         **
!**                                                         **
!*************************************************************
!*************************************************************

9000 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent environment file : '//envfich
ERRVAL=-1
goto 9999

9005 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid environment file : '//envfich
ERRVAL=-2
goto 9999

9010 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid name for the environment file'
ERRVAL=-3
goto 9999

9015 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid name for the list file'
ERRVAL=-4
goto 9999

9020 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid name for the scene file'
ERRVAL=-5
goto 9999

9025 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid name for the physical characteristic file'
ERRVAL=-6
goto 9999

9030 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid name for the geoarea file'
ERRVAL=-7
goto 9999

9035 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid name for the geopix file'
ERRVAL=-8
goto 9999

9040 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid name for the geovisil file'
ERRVAL=-9
goto 9999

9045 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid name for the spectral characteristic file'
ERRVAL=-10
goto 9999

9050 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent geovisil file : '//geovisil
ERRVAL=-11
goto 9999

9055 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent geopix file : '//geopix
ERRVAL=-12
goto 9999

9060 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid geovisil file : '//geovisil
ERRVAL=-13
goto 9999

9065 ERRKEY = 'SPECTRIMAG_MEM_ERR'
err_msg='Not enough memory ressources to run the program; allocation '//err_code_mem
ERRVAL=-14
goto 9999

9070 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid geopix file : '//geopix
ERRVAL=-15
goto 9999

9075 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent scene file : '//scenarea
ERRVAL=-16
goto 9999

9080 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent physical characteristic file : '//charactarea
ERRVAL=-17
goto 9999

9085 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent geoarea file : '//geoarea
ERRVAL=-18
goto 9999

9090 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent list file : '//listarea
ERRVAL=-19
goto 9999

9095 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent spectral characteristic file : '//charaspectr
ERRVAL=-20
goto 9999

9100 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid unit for the input optical files'
ERRVAL=-21
goto 9999

9105 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid unit for the output cube file'
ERRVAL=-22
goto 9999

9110 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid spectral characteristic file : '//charaspectr
ERRVAL=-23
goto 9999

9111 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid image option : '//opt3
ERRVAL=-23
goto 9999

9112 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid correction option : '//opt4
ERRVAL=-23
goto 9999

9113 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid interpolation option : '//opt5
ERRVAL=-23
goto 9999

9115 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid scene file : '//scenarea
ERRVAL=-24
goto 9999

9116 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid area option : '//opt1
ERRVAL=-24
goto 9999

9117 ERRKEY = 'SPECTRIMAG_I/O_ERR'
write(err_msg,FMT='("Too many compounds in layer ",I1," of area ",I1000)')&
    &m,nbta
ERRVAL=-24
goto 9999

9118 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid texture : '//tex(m)
ERRVAL=-24
goto 9999

9119 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid optical depth option : '//opt6
ERRVAL=-24
goto 9999

9120 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent optical file : '//name
ERRVAL=-25
goto 9999

9125 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid optical file : '//name
ERRVAL=-26
goto 9999

9130 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid physical characteristic file : '//charactarea
ERRVAL=-27
goto 9999

9135 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid list file : '//listarea
ERRVAL=-28
goto 9999

9140 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid geoarea file : '//geoarea
ERRVAL=-29
goto 9999

9145 ERRKEY = 'SPECTRIMAG_RANGE_ERR'
err_msg='User spectral domain not represented in the optical files'
ERRVAL=-30
goto 9999

9150 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent deresolution file : '//derefich
ERRVAL=-31
goto 9999

9155 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid deresolution file : '//derefich
ERRVAL=-32
goto 9999

9160 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Error writing the backplane information'
ERRVAL=-33
goto 9999

9165 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Error writing the BAND_BIN group'
ERRVAL=-34
goto 9999

9170 ERRKEY = 'SPECTRIMAG_DERE_ERR'
err_msg='Spectral convolution error : '//err_msg
ERRVAL=-35
goto 9999

9175 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Scratch files already exist'
ERRVAL=-36
goto 9999

9176 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Read or write error on the scratch files'
ERRVAL=-37
goto 9999

9177 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Non existent atmosphere file : '//atmtau
ERRVAL=-38
goto 9999

9178 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Invalid atmosphere file : '//atmtau
ERRVAL=-39
goto 9999

9179 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Raw file already exists'
ERRVAL=-40
goto 9999

9180 ERRKEY = 'SPECTRIMAG_I/O_ERR'
err_msg='Raw file cannot be created if number of blocks exceeds 1'
ERRVAL=-40
goto 9999

!****************************************************************************
! Write an error message and return to DOCUBE to inform of error
!***************************************************************************
9999 continue

WRITE(UNIT=6,FMT='(A19," : ",A255)')ERRKEY,err_msg

contains

function fctP(g,alpha)
 implicit none
 real(kind=8):: fctP, g, alpha
 fctP=(1-g**2)/(1+g**2-2.*g*cos(alpha))**(1.5)
 return
end function fctP

end program







