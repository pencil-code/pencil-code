;-------------------------------------------------------------
;+
; NAME:
;       FACTOR
; PURPOSE:
;       Find prime factors of a given number.
; CATEGORY:
; CALLING SEQUENCE:
;       factor, x, p, n
; INPUTS:
;       x = Number to factor (>1).       in
; KEYWORD PARAMETERS:
;       Keywords:
;         /QUIET  means do not print factors.
;         /DEBUG  Means list steps as they happen.
;         /TRY    Go beyond 20000 primes.
; OUTPUTS:
;       p = Array of prime numbers.      out
;       n = Count of each element of p.  out
; COMMON BLOCKS:
; NOTES:
;       Note: see also prime, numfactors, print_fact.
; MODIFICATION HISTORY:
;       R. Sterner.  4 Oct, 1988.
;       RES 25 Oct, 1990 --- converted to IDL V2.
;       R. Sterner, 1999 Jun 30 --- Improved (faster, bigger).
;       R. Sterner, 1999 Jul  7 --- Bigger values (used unsigned).
;       R. Sterner, 1999 Jul  9 --- Tried to make backward compatable.
;       R. Sterner, 2000 Jan 06 --- Fixed to ignore non-positive numbers.
;       Johns Hopkins University Applied Physics Laboratory.
;
; Copyright (C) 1988, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
; NAME:
;       SPC
; PURPOSE:
;       Return a string with the specified number of spaces (or other char).
; CATEGORY:
; CALLING SEQUENCE:
;       s = spc(n, [text])
; INPUTS:
;       n = number of spaces (= string length).   in 
;        text = optional text string.              in
;          # spaces returned is n-strlen(strtrim(text,2))
; KEYWORD PARAMETERS:
;       Keywords:
;         CHARACTER=ch  Character other than a space.
;           Ex: CHAR='-'.
;         /NOTRIM means do not do a strtrim on text.
; OUTPUTS:
;       s = resulting string.                     out
; COMMON BLOCKS:
; NOTES:
;       Note: Number of requested spaces is reduced by the
;         length of given string.  Useful for text formatting.
; MODIFICATION HISTORY:
;       Written by R. Sterner, 16 Dec, 1984.
;       RES --- rewritten 14 Jan, 1986.
;       R. Sterner, 27 Jun, 1990 --- added text.
;       R. Sterner, 1994 Sep  7 --- Allowed text arrays.
;       R. Sterner, 1999 Jul  2 --- Added /NOTRIM keyword.
;       Johns Hopkins University Applied Physics Laboratory.
;
; Copyright (C) 1984, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-------------------------------------------------------------
 
	function spc,n, text, character=char, notrim=notrim, help=hlp
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Return a string with the specified number of spaces (or '+$
	    'other char).' 
	  print,' s = spc(n, [text])' 
	  print, '  n = number of spaces (= string length).   in '
	  print,'   text = optional text string.              in'
	  print,'     # spaces returned is n-strlen(strtrim(text,2))'
	  print,'   s = resulting string.                     out' 
	  print,' Keywords:'
	  print,'   CHARACTER=ch  Character other than a space.'
	  print,"     Ex: CHAR='-'."
	  print,'   /NOTRIM means do not do a strtrim on text.'
	  print,' Note: Number of requested spaces is reduced by the'
	  print,'   length of given string.  Useful for text formatting.'
	  return, -1
	endif
 
	if n_params(0) eq 1 then begin
	  n2 = n
	endif else begin
	  if keyword_set(notrim) then $
	    ntxt=strlen(text) else ntxt=strlen(strtrim(text,2))
;	  n2 = n - strlen(strtrim(text,2))
	  n2 = n - ntxt
	endelse
 
	ascii = 32B
	if n_elements(char) ne 0 then ascii = (byte(char))[0]
 
	num = n_elements(n2)
	out = strarr(num)
	for i = 0, num-1 do begin
	  if n2[i] le 0 then out[i] = '' else $
	    out[i] = string(bytarr(n2[i]) + ascii)
	endfor
 
	if n_elements(out) eq 1 then out=out[0]
	return, out
 
	end


;-------------------------------------------------------------
; NAME:
;       PRINT_FACT
; PURPOSE:
;       Print prime factors found by the factor routine.
; CATEGORY:
; CALLING SEQUENCE:
;       print_fact, p, n
; INPUTS:
;       p = prime factors.          in
;       n = number of each factor.  in
; KEYWORD PARAMETERS:
; OUTPUTS:
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       R. Sterner  4 Oct, 1988.
;       RES 25 Oct, 1990 --- converted to IDL V2.
;       R. Sterner, 26 Feb, 1991 --- Renamed from print_factors.pro
;       R. Sterner, 1999 Jun 30 --- Better output format.
;       R. Sterner, 1999 Jul  7 --- Bigger values (used unsigned).
;       R. Sterner, 1999 Jul  9 --- Made backward compatable.
;
; Copyright (C) 1988, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-------------------------------------------------------------
 
	pro print_fact, p, n, help=hlp
 
	if (n_params(0) lt 2) or keyword_set(hlp) then begin
	  print,' Print prime factors found by the factor routine.'
	  print,' print_fact, p, n'
	  print,'   p = prime factors.          in'
	  print,'   n = number of each factor.  in'
	  return
	endif
 
	;-------  Drop unused primes  ---------------
	w = where(n gt 0)	; Find only primes used.
	p2 = p[w]
	n2 = n[w]
 
	;-------  Use largest available integer type  --------------
	flag = !version.release ge 5.2
	if flag eq 1 then begin
	  err=execute('t=1ULL')		; Use 64 bit int (hide from old IDL).
	endif else begin
	  t = 1L			; Use long int (best available in old).
	endelse
 
	;-------  Compute number from it's prime factors.  ----------
	for i = 0, n_elements(p2)-1 do t = t * p2[i]^n2[i]
 
	;-------  Prepare output  -----------------------
	a = strtrim(t,2)+' = '			; Start factors string.
	b = ''					; Start exponents string.
	last = n_elements(p2)-1			; Last factors index.
	for i=0, last do begin
	  a = a + strtrim(p2[i],2)		; Insert next factor.
	  lena = strlen(a)			; Length of factor string.
	  nxtb = strtrim(n2[i],2)		; Next exponent.
	  if nxtb eq '1' then nxtb=' '		; Weed out 1s.
	  b = b+spc(lena,b,/notrim)+nxtb	; Insert next exponent.
	  if i ne last then a=a+' x '		; Not last, add x.
	endfor
 
	;------  Print exponents and factors  -----------
	print,' '
	print,' '+b
	print,' '+a
 
	return
	end


 
	pro factor, x, p, n, quiet=quiet, debug=debug, try=try, help=hlp
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Find prime factors of a given number.'
	  print,' factor, x, p, n'
	  print,'   x = Number to factor (>1).       in'
	  print,'   p = Array of prime numbers.      out'
	  print,'   n = Count of each element of p.  out'
	  print,' Keywords:'
	  print,'   /QUIET  means do not print factors.'
	  print,'   /DEBUG  Means list steps as they happen.'
	  print,'   /TRY    Go beyond 20000 primes.'
	  print,' Note: see also prime, numfactors, print_fact.'
	  return
	endif
 
	if x le 0 then return
 
	flag = !version.release ge 5.2
 
	s = sqrt(x)			; Only need primes up to sqrt(x).
	g = long(50 + 0.13457*s)	; Upper limit of # primes up to s.
	np = 50				; Start with np (50) primes.
	p = prime(np)			; Find np primes.
	n = intarr(n_elements(p))	; Divisor count.
 
	if flag eq 1 then $		; Working number.
	  err=execute('t=ulong64(x)') $	; Use best integer available.
	  else t=long(x)		; Best pre-5.2 integer.
	i = 0L				; Index of test prime.
 
loop:	pt = p[i]			; Pull test prime.
	if keyword_set(debug) then $
	  print,' Trying '+strtrim(pt,2)+' into '+strtrim(t,2)
	if flag eq 1 then $
	  err=execute('t2=ulong64(t/pt)') $
	  else t2=long(t/pt)
	if t eq t2*pt then begin	; Check if it divides.
	  if keyword_set(debug) then $
	    print,'   Was a factor.  Now do '+strtrim(t2,2)
	  n[i] = n[i] + 1		; Yes, count it.
	  t = t2			; Result after division.
	  if t2 eq 1 then goto, done	; Check if done.
	  goto, loop			; Continue.
	endif else begin
	  i = i + 1			; Try next prime.
	  if i ge np then begin
	    s = sqrt(t)			; Only need primes up to sqrt(x).
	    g = long(50 + 0.13457*s)	; Upper limit of # primes up to s.
	    if g le np then goto, last	; Must be done.
	    np = (np+50)<g		; Want 50 more primes.
	    if (np gt 20000) and (not keyword_set(try)) then begin
	      print,' Too hard.  Tried '+strtrim(np-50,2)+' primes.'
	      print,' Trying to crack '+strtrim(t,2)
	      print,' To go farther use keyword /TRY.'
	      return
	    endif
	    if keyword_set(debug) then $
	      print,' Finding more primes: '+strtrim(np,2)+ $
	      '.  Max needed = '+strtrim(g,2)
	    p = prime(np)		; Find primes.
	    n = [n,intarr(50)]		; Make room for more factors.
	  endif
	  if i ge g then goto, last	; Nothing up to sqrt works.
	  goto, loop			; Continue.
	endelse
 
last:	p = [p,t]			; Residue was > sqrt, must be prime.
	n = [n,1]			; Must occur only once. (else < sqrt).
 
done:	w = where(n gt 0)
	n = n[w]			; Trim excess off tables.
	p = p[w]
 
	if not keyword_set(quiet) then print_fact, p, n
 
	return
	end

