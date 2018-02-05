pro _cache_index_structures

; don't waste time reading in the same index tables repeatedly ...

  COMMON astrom_index_tables, atlas, index
  if n_elements(index) EQ 0 then begin
      ; bad to have hardcoded file paths, could use env variable in future
      atlas = mrdfits('../etc/astrom-atlas.fits', 1)
      index = mrdfits('../etc/tr_neo2_index.fits', 1)
  endif

end

function nominal_astrometry, coadd_id

  _cache_index_structures
  COMMON astrom_index_tables, atlas, index

  ; this is slow, could cache sorted version of tables and binary search
  w = (where(atlas.coadd_id EQ coadd_id, nw))[0]

  if (nw NE 1) then stop ; sanity check

  return, atlas[w]
end

function recalib_astrometry, coadd_id, epoch, band

  _cache_index_structures
  COMMON astrom_index_tables, atlas, index

  ; this is slow, could cache sorted version of tables and binary search
  w = (where((index.coadd_id EQ coadd_id) AND (index.epoch EQ epoch) AND $
             (index.band EQ band), nw))[0]

  if (nw NE 1) then stop

  return, index[w]
end

function compute_pixel_offsets, coadd_id, epoch, band

  astr_orig = nominal_astrometry(coadd_id) 
  astr_new = recalib_astrometry(coadd_id, epoch, band)

; could also make these once and then used cached version as an optimization
  xbox = lindgen(astr_orig.naxis[0], astr_orig.naxis[1]) MOD astr_orig.naxis[0]
  ybox = lindgen(astr_orig.naxis[0], astr_orig.naxis[1]) / astr_orig.naxis[0]

; ** this is the actual computation of interest **
  xy2ad, xbox, ybox, astr_orig, abox, dbox
  ad2xy, abox, dbox, astr_new, xx, yy

; could use median instead i suppose ...
  dx = mean(xx-double(xbox))
  dy = mean(yy-double(ybox))

  return, [dx, dy]
end

function apply_sinc_shift, im, coadd_id, epoch, band

  if n_elements(coadd_id) NE 1 then stop
  if n_elements(epoch) NE 1 then stop

  sz = size(im, /dim)

  if n_elements(sz) NE 2 then stop
  if (sz[0] NE 2048) OR (sz[1] NE 2048) then stop

  offs = compute_pixel_offsets(coadd_id, epoch, band)

; shift the image so as to make it roughly registered properly according
; to the nominal WCS
  ims = sshift2d(im, -1.0d*offs)

  return, ims
end

pro _example, diff_raw, diff_shifted

  coadd_id = '2246p106'
  band = 1

  cube_raw = fltarr(2048, 2048, 2)
  cube_shifted = fltarr(2048, 2048, 2)
  for epoch=5, 6 do begin
      fname = $
        '/global/projecta/projectdirs/cosmo/work/wise/outputs/merge/neo2/e' + $
         string(epoch, format='(I03)') + '/' + strmid(coadd_id, 0, 3) + '/' + $
         coadd_id + '/unwise-' + coadd_id + '-w1-img-u.fits'
      im = readfits(fname)
      cube_raw[*, *, epoch-5] = im
      cube_shifted[*, *, epoch-5] = apply_sinc_shift(im, coadd_id, epoch, band)
  endfor

  diff_raw = cube_raw[*,*,1] - cube_raw[*, *, 0]
  diff_shifted = cube_shifted[*,*,1] - cube_shifted[*, *, 0]

end
