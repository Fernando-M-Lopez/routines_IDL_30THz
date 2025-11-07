pro correction_jitter_30THz, files, sector=sector, Iref=Iref, gif=gif, sav=sav, wavelet=wavelet, factor = factor, outdir = outdir

	; ****************************************Information ********************************************************** 
	; The routine realize the correction of the jitter in the 30THz images for the AR30T and BR30T telescopes
	; The correction is performed by applying a Cross-crorrelation between the image of interest and a reference one. 
	; it is necessary that the images by previously calibrated by Flat field and centered. 
	; The application of a wavelet Atrous can be used to eliminate the noise in the images
	; 
	; Procedure: 
	; The routine will require to select a ROI with a clear structure (e.g. Active Region) to perform the correlation
	; The aligned images are save as level2.0 fits files 
	; 
	; call: 
	; correction_jitter_30THz, files
	; 
	; input: 
	; files       = 'array of strings with the list of images to correct' 
	; ROI 	= [x0, y0, nx, ny] the left corner and size of the ROI to analize - If not indicated the program will require to select a box   
	;	
	; output: 
	; the images aligned by Cross-Correlation
	; 
	; Optional: 
	; Iref=X        = The number (integer) of the reference image in the list of files (if ommited the Iref=0 is assumed) 
	; ,/wavelet 	  = to eliminate the noise in the images using discrete Atrous wavelet deconposition 
	; ,/gif         = this keyword will save the ROI of all images as gif files (not done yet) 
	; ,/sav         = use this keyword to save the ROI sector for all images in a .sav file (not done yet)
	; factor = X     = To correct subpixels scale, then rebin the ROI to X times the original size
	; 
	; History:
	; Written by Fernando M. Lopez (CRAAM-Mackenzie) --- November 2019
	;*****************************************************************************************************************
	loadct,0
		if n_elements(outdir) eq 0 then begin 
			pos_barra=strpos(files[0],'/',/reverse_search)
			outdir=strmid(files[0],0,pos_barra+1) + 'level20/'
			file_mkdir, outdir
		endif
		
	; read the images

	if size(files, /tname) ne 'STRING' then message, '***** ERROR: THE ROUTINE REQUIERE FILE LIST AS STRING *******'
	nimages = n_elements(files)

	;*****************************************************************************************************************
	; read the first image and select the ROI

	Image0 = readfits(files[0], header0)
	ss = size(Image0)
	ximage = ss[1] & yimage = ss[2]
	
	if n_elements(sector) ne 4 then begin 
		message, '******* Must select a sector = [x0,y0,nx,ny] ********', /continue
		window,0, xs=ximage/2, ys=yimage/2
		tvscl, sqrt(rebin(Image0,ximage/2,yimage/2))
		box_cursor, x0,y0,nx,ny
		sector = fltarr(4)
		sector[0]=2*x0 & sector[1]=2*y0 & sector[2]=nx*2 & sector[3]=ny*2
	endif 

	x0=sector[0] & y0=sector[1] & nx=sector[2] & ny=sector[3]


	dummy23=Image0

	dummy23[x0-1:x0,y0:y0+ny-1]=10000.
	dummy23[x0:x0+nx-1,y0-1:y0]=10000.
	dummy23[x0+nx-1:x0+nx,y0:y0+ny-1]=10000.
	dummy23[x0:x0+nx-1,y0+ny-1:y0+ny]=10000. 

	;*****************************************************************************************************************
	; Show the ROI in the original image to check it
	window,1, xs=700, ys=700
	wset,1 & tvscl, sqrt(congrid(dummy23,700,700,/cubic))

	;*****************************************************************************************************************
	; Read all the images and start the cross-correlation

	; cube of ROI with the aligned images to save!!!
	ROI = fltarr(nx, ny, nimages) 
	ROI[*,*,0] = Image0[x0:x0+nx-1,y0:y0+ny-1]

	headers = strarr(n_elements(header0)+1, nimages)
	
	valor_rms = fltarr(nimages, 3) ; valor_rms(*,0) = rms original, valor_rms(*,1) = rms metodo 1, valor_rms(*,2) = rms metodo 2   
	
	If ~keyword_set(factor) then factor = 1.
	If ~keyword_set(Iref) then Iref=0

	for i = 0l, nimages-1 do begin
		; reference image	
		Image0_orig = readfits(files[Iref], header0)
		if keyword_set(wavelet) then Image0_orig = filtered_atrous(Image0_orig)
		region0_orig = Image0_orig[x0:x0+nx-1,y0:y0+ny-1] 

		Image0_orig_atrous = Image0_orig
		Image0_atrous = Image0_orig_atrous
		
		Image0_atrous = Image0_atrous - min(Image0_atrous)
		Image0_atrous = Image0_atrous / max(Image0_atrous)
		region0_atrous = Image0_atrous[x0:x0+nx-1,y0:y0+ny-1]
		region0_atrous = Image0_atrous[x0:x0+nx-1,y0:y0+ny-1] 

		; Image to align 			
		ImageI_orig = readfits(files[i],headeri)
		if keyword_set(wavelet) then ImageI_orig = filtered_atrous(ImageI_orig)

		ImageI_orig_shifted = ImageI_orig
		regionI_orig = ImageI_orig[x0:x0+nx-1,y0:y0+ny-1]	 

		ImageI_orig_atrous = ImageI_orig	
		
		ImageI_atrous = ImageI_orig_atrous

		Imagei_atrous = Imagei_atrous - min(Imagei_atrous)
		Imagei_atrous = Imagei_atrous / max(Imagei_atrous)
		regionI_atrous = ImageI_atrous[x0:x0+nx-1,y0:y0+ny-1]	 

	;the regions to correlate are: 
	; region0_atrous & regionI_atrous



	;################################################# METHOD 1 FOR CORRELATION ##########################################################

	lag1 = 10 ; the value of pixels to +- shift [lagx,lagy]

	; correct integer pixels in the image
		Result11 = CROSS_CORR2_CAR(region0_atrous,regionI_atrous,lag1, /REPORT)
		ImageI_orig_shifted = shift(ImageI_orig_shifted, ROUND(-1.* Result11(0)), ROUND(-1.* Result11(1)))  
		; region shifted by the 1st approach:
		regionI_orig_shifted11 = ImageI_orig_shifted[x0:x0+nx-1,y0:y0+ny-1]	
		ImageI_corr1_atrous = shift(ImageI_atrous, ROUND(-1.* Result11(0)), ROUND(-1.* Result11(1)))

	;To correct subpixels scale, then rebin the ROI to "factor" times the original size
		regionI_atrous = ImageI_corr1_atrous[x0:x0+nx-1,y0:y0+ny-1]	
		result12 = ALIGN_CROSS_CAR(region0_atrous,regionI_atrous)
		ImageI_orig_shifted_big = REBIN(ImageI_orig_shifted,ximage*factor,yimage*factor)	
		ImageI_orig_shifted_big = shift(ImageI_orig_shifted_big, ROUND(1.* Result12(0)*factor), ROUND(1.* Result12(1)*factor))

		ImageI_orig_shifted12 = rebin(ImageI_orig_shifted_big,ximage,yimage)
		
		;  region shifted by the 2nd approach:
		regionI_orig_shifted12 = ImageI_orig_shifted12[x0:x0+nx-1,y0:y0+ny-1]

		if rms(region0_orig - regionI_orig) gt 10. then print, ' TO HIGH RMS = ',rms(region0_orig - regionI_orig),' SKIPPING CORRELATION'
		if rms(region0_orig - regionI_orig) gt 10. then continue
		

		; Calculates the Root-mean-square for the (region0 - regionI_corr)
		print, '****************************************************************************************************'
		print, ' Correlating image number: ', i, ' of a total of: ', nimages
		print, '****************************************************************************************************'
		print, ' Writting:                                                    RMS  '
		print, 'Root-mean-square of the Original data rms(Ii - Iref)= ', rms(region0_orig - regionI_orig)
		print, '****************************************************************************************************'

		; ----------------------------------------- for method 1 ---------------------------------------------------
		; Calculates the residuals for the (region0 - regionI_corr)
		print, ' For method 1 '
		print, 'Root-mean-square of the Correlation for 1st aproach = ', rms(region0_orig - regionI_orig_shifted11)
		print, 'Root-mean-square of the Correlation for 2nd aproach = ', rms(region0_orig - regionI_orig_shifted12)
		print, 'shift final del metodo 1: x=', result11[0],' + ',result12[0], ', y=',result11[1],' + ',result12[1]
		print, '****************************************************************************************************'
		; -----------------------------------------------------------------------------------------------------------



	;################################################# METHOD 2 FOR CORRELATION ##########################################################

	; CORREL_IMAGES:
	; MAGNIFICATION = option causes computation of high resolution correlation
	;                      of magnified images, thus much slower.
	;                      Shifting distance is automatically = 2 + Magnification,
	;                      and optimal pixel offset should be known and specified.
	;                      Optimal offset can then be found to fractional pixels
	;                      using CorrMat_Analyze( correl_images( ) ).

		; reinitializate the values of the Images


		; the original image 0: Image0_orig
		; the original image I: ImageI_orig

	 ImageI_atrous = ImageI_orig_atrous
	 Image0_atrous = Image0_orig_atrous

		ImageI_atrous = ImageI_atrous - min(ImageI_atrous)
		ImageI_atrous = ImageI_atrous / max(ImageI_atrous)
		Image0_atrous = Image0_atrous - min(Image0_atrous)
		Image0_atrous = Image0_atrous / max(Image0_atrous)

		regionI_atrous = ImageI_atrous[x0:x0+nx-1,y0:y0+ny-1]	
		region0_atrous = Image0_atrous[x0:x0+nx-1,y0:y0+ny-1]

		; shift the image I pixels dimension
		result21 = correl_images(region0_atrous, regionI_atrous)
		corrmat_analyze, result21, xoffset21, yoffset21
		ImageI_corr21_atrous = shift(ImageI_atrous, round(1.*xoffset21), round(1.*yoffset21)) 
		regionI_corr21_atrous = ImageI_corr21_atrous[x0:x0+nx-1,y0:y0+ny-1]
		
		; the original Image I shifted: 
		ImageI_shifted21 = shift(ImageI_orig, round(1.*xoffset21), round(1.*yoffset21))
		regionI_corr21 = ImageI_shifted21[x0:x0+nx-1,y0:y0+ny-1]

		; shift the image I sub-pixels dimension
		result22 = correl_images(region0_atrous, regionI_corr21_atrous, magnification=factor)
		corrmat_analyze, result22, xoffset22, yoffset22,magnification=factor
		ImageI_shifted22_big = rebin(ImageI_shifted21,ximage*factor,yimage*factor)
		ImageI_shifted22_big = shift(ImageI_shifted22_big, round(1.*xoffset22*factor), round(1.*yoffset22*factor))

		ImageI_shifted22 = rebin(ImageI_shifted22_big, ximage, yimage)

		regionI_corr22 = ImageI_shifted22[x0:x0+nx-1,y0:y0+ny-1]

		; ----------------------------------------- for method 2 ---------------------------------------------------
		; Calculates the Root-mean-square for the ROI where (region0 - regionI_corr22)
		print, ' For method 2 '
		print, 'Root-mean-square of the Correlation for 1st aproach = ', rms(region0_orig - regionI_corr21)
		print, 'Root-mean-square of the Correlation for 2nd aproach = ', rms(region0_orig - regionI_corr22)
		print, 'shift final del metodo 2: x=',xoffset21 ,' + ',xoffset22, ', y=',yoffset21,' + ',yoffset22
		print, '****************************************************************************************************'
		; -----------------------------------------------------------------------------------------------------------

			valor_rms[i,0] = rms(region0_orig - regionI_orig) & valor_rms[i,1] = rms(region0_orig - regionI_orig_shifted12) & $
						valor_rms[i,2] = rms(region0_orig - regionI_corr22)

	 	if rms(region0_orig - regionI_orig) lt rms(region0_orig - regionI_orig_shifted12) and $ 
						rms(region0_orig - regionI_orig) lt rms(region0_orig - regionI_corr22) then begin

		message, '#### THE CROSS-CORRELATION APPROACH FAILED: MUST USE AN ACTIVE REGION TO ALIGN THE IMAGES #####',/CONTINUE
				ImageI_shifted_save = Image0_orig
		endif

		if rms(region0_orig - regionI_orig_shifted12) ge rms(region0_orig - regionI_corr22) then ImageI_shifted_save = ImageI_shifted22 $
			else  ImageI_shifted_save = ImageI_orig_shifted12

		; ********************************************************************************************************************************
		; save the full aligned images as a level 2.0 data

		filename = ''
		filename = outdir + strmid(sxpar(headeri, 'FILENAME'),0,35)+'20.fits'
		sxaddpar,headeri,'FILENAME',strmid(sxpar(headeri, 'FILENAME'),0,35)+'20.fits'
		sxaddpar,headeri,'LVL_NUM', '2.0'		
		sxaddpar, headeri, 'HISTORY', 'Image aligned using active region by Cross-Correlation technique'
		writefits, filename, ImageI_shifted_save, headeri

		; ********************************************************************************************************************************

		; ROI is the region of interest to save!!

		
		ROI[*,*,i] = ImageI_shifted_save[x0:x0+nx-1,y0:y0+ny-1] 
	 
		headers[*,i] = headeri[0:n_elements(header0)]
	
	endfor

	; save the information in a .sav file

	savename = strmid(sxpar(header0,'TELESCOP'),0,5)+'-'+ strmid(sxpar(header0,'DATE_OBS'),0,10) + '-ROI-aligned-noise_wavelet.sav'
	
	if keyword_set(sav) then save, ROI, headers, sector, valor_rms, description = '*** ROI is the aligned ROI of each image; headers are the headers of the images; sector = is the [x0,y0, Nx, Ny] location of the ROI in the 800x800 30THz images*** ', filename=savename1, /verbose 


stop
end


	; ############################################################################################################################################
	function filtered_atrous, Image
	; This routine apply a atrus wavelet to enhance the structure of the active regions to improve the cross-correlation technique
	; The routine eliminate noise in all the scales, and eliminate from the image the smoothed scale!!!!
	; Description of the atrous wavelet in Starck & Murtagh, "Astronomical Image anda Data Analysis" - Second Edition. 

	sz = size(image)
	filter = [1/16., 1/4., 3/8.,1/4., 1/16.]
	;nscales = floor(alog((sz[1] < sz[2])/n_elements(filter))/alog(2))

	nscales=5

	npix = n_elements(Image) ; number of pixels of the image
	sigma = stdev(Image)

	fmat = filter#transpose(filter)

	decomp = fltarr(sz[1], sz[2], nscales+1)
	 im = image
	for j = 0, nscales-1 do begin

		smooth = convolve(im, fmat)
		decomp[*,*,nscales-j] = im - smooth
		im = smooth
		; Generate new filter
		  newfilter = fltarr(2*n_elements(filter)-1) 
		  newfilter[2*findgen(n_elements(filter))] = filter
	; filter is padded with zeros between the images
		  fmat = newfilter#transpose(newfilter)
		  filter = newfilter
	endfor

	decomp[*,*,0] = smooth

	sigma = fltarr(nscales+1)

	for i=0,nscales -1 do begin
		sigma[nscales-i] = stdev(decomp[*,*,nscales-i])
		thredshold = sqrt(2.*alog(npix)) * sigma  
; ;		the threshold to eliminate Gaussian noise in the wavelet scales 
; (Stark & Murtagh)
		result = where(decomp[*,*,nscales-i] lt thredshold[nscales-i],count)
		if count gt 0 then decomp[result,nscales-i] = 0.
	endfor

	image_filtered = decomp[*,*,2] + decomp[*,*,3] + decomp[*,*,4] + decomp[*,*,5] + decomp[*,*,0] + decomp[*,*,1]
	return, image_filtered
	end




