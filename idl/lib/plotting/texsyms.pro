;;;;;;;;;;;;;;;;;;;;;;;      -*- coding:latin1 -*-  vim:set encoding=latin1:
;;;   texsyms.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   27-Mar-2002
;;;
;;;  Description:
;;;   Symbolic names for symbols that work for X11 as well as for
;;;   PostScript. Supposedly more convenient than textoidl with its
;;;   function calls. The nomenclature for the symbols is <TeX-name>
;;;   for lowercase letters and <TeX-name>_c fur capital letters.
;;;  Usage:
;;;   s = texsyms()
;;;   plot, .., TITLE=s.alpha+'!6-'+s.omega_c+'!6-dynamo!X'

function texsyms, DEVICE=device

  default, device, !d.name

  if ((device eq 'PS') and (!p.font eq 0)) then begin
    ; lowercase lsymbols
    sym = { $
              alpha      : '!7a!X' , $
              beta       : '!7b!X' , $
              gamma      : '!7g!X' , $
              delta      : '!7d!X' , $
              epsilon    : '!7e!X' , $
              zeta       : '!7z!X' , $
              eta        : '!7h!X' , $
              theta      : '!7q!X' , $
              iota       : '!7i!X' , $
              kappa      : '!7k!X' , $
              lambda     : '!7l!X' , $
              mu         : '!7m!X' , $
              nu         : '!7n!X' , $
              xi         : '!7x!X' , $
              omicron    : '!7o!X' , $
              pi         : '!7p!X' , $
              rho        : '!7r!X' , $
              sigma      : '!7s!X' , $
              tau        : '!7t!X' , $
              ypsilon    : '!7u!X' , $
              phi        : '!7f!X' , $
              chi        : '!7c!X' , $
              psi        : '!7y!X' , $
              omega      : '!7w!X' , $
              varepsilon : '!7e!X' , $
              vartheta   : '!7J!X' , $
              varpi      : '!7v!X' , $
              varrho     : '!7r!X' , $
              varphi     : '!7j!X' , $
    ; uppercase symbols
              alpha_c     : '!7A!X' , $
              beta_c      : '!7B!X' , $
              gamma_c     : '!7G!X' , $
              delta_c     : '!7D!X' , $
              epsilon_c   : '!7E!X' , $
              zeta_c      : '!7Z!X' , $
              eta_c       : '!7H!X' , $
              theta_c     : '!7Q!X' , $
              iota_c      : '!7I!X' , $
              kappa_c     : '!7K!X' , $
              lambda_c    : '!7L!X' , $
              mu_c        : '!7M!X' , $
              nu_c        : '!7N!X' , $
              xi_c        : '!7X!X' , $
              omicron_c   : '!7O!X' , $
              pi_c        : '!7P!X' , $
              rho_c       : '!7R!X' , $
              sigma_c     : '!7S!X' , $
              tau_c       : '!7T!X' , $
              ypsilon_c   : '!7U!X' , $
              phi_c       : '!7F!X' , $
              chi_c       : '!7C!X' , $
              psi_c       : '!7Y!X' , $
              omega_c     : '!7W!X' , $
    ; other symbols
              cdot        : '!6·!X' , $
              mid         : '!6|!X' , $
              Im          : '!7Á!X' , $
              Re          : '!7Â!X' , $
              sqrt        : '!7Ö!X' , $
              times       : '!7´!X' , $
              leq         : '!6!S<!R_!X' , $ ; `le' is reserved..
              geq         : '!6!S>!R_!X' , $ ; `ge' is reserved..
              int         : '!Mò!X' , $
              partial     : '!M¶!X' , $
              nabla       : '!MÑ!X' , $
              uparrow     : '!7­!X', $
              downarrow   : '!7¯!X', $
              leftarrow   : '!7¬!X', $
              rightarrow  : '!7®!X', $
              infty       : '!7¥!X',  $
    ; harder stuff like tilde on top of symbol, etc
    ; Usage: PLOT, [...], TITLE=sym.tilde+'!8Q!X' to get \tilde{Q}, etc
              tilde       : '!S!N!A!3~!R!X', $
              bar         : '!S!N!A!3-!R!X' $
  }
  endif else begin ; Hershey fonts
    ; lowercase lsymbols
    sym = { $
              alpha      : '!7a!X' , $
              beta       : '!7b!X' , $
              gamma      : '!7c!X' , $
              delta      : '!7d!X' , $
              epsilon    : '!7e!X' , $
              zeta       : '!7f!X' , $
              eta        : '!7g!X' , $
              theta      : '!7h!X' , $
              iota       : '!7i!X' , $
              kappa      : '!7j!X' , $
              lambda     : '!7k!X' , $
              mu         : '!7l!X' , $
              nu         : '!7m!X' , $
              xi         : '!7n!X' , $
              omicron    : '!7o!X' , $
              pi         : '!7p!X' , $
              rho        : '!7q!X' , $
              sigma      : '!7r!X' , $
              tau        : '!7s!X' , $
              ypsilon    : '!7t!X' , $
              phi        : '!7u!X' , $
              chi        : '!7v!X' , $
              psi        : '!7w!X' , $
              omega      : '!7x!X' , $
              varepsilon : '!7e!X' , $
              vartheta   : '!7h!X' , $
              varpi      : '!S!7x!R!A_!N!X' , $
              varrho     : '!7q!X' , $
              varphi     : '!7u!X' , $
    ; uppercase symbols
              alpha_c     : '!7A!X' , $
              beta_c      : '!7B!X' , $
              gamma_c     : '!7C!X' , $
              delta_c     : '!7D!X' , $
              epsilon_c   : '!7E!X' , $
              zeta_c      : '!7F!X' , $
              eta_c       : '!7G!X' , $
              theta_c     : '!7H!X' , $
              iota_c      : '!7I!X' , $
              kappa_c     : '!7J!X' , $
              lambda_c    : '!7K!X' , $
              mu_c        : '!7L!X' , $
              nu_c        : '!7M!X' , $
              xi_c        : '!7N!X' , $
              omicron_c   : '!7O!X' , $
              pi_c        : '!7P!X' , $
              rho_c       : '!7Q!X' , $
              sigma_c     : '!7R!X' , $
              tau_c       : '!7S!X' , $
              ypsilon_c   : '!7T!X' , $
              phi_c       : '!7U!X' , $
              chi_c       : '!7V!X' , $
              psi_c       : '!7W!X' , $
              omega_c     : '!7X!X' , $
    ; other symbols
              cdot        : '!3·!X',  $
              mid         : '!3|!X' , $
              Im          : '!15I!X', $
              Re          : '!15R!X', $
              sqrt        : '!9r!X',  $
              times       : '!9X!X',  $
              odot        : '!9n!X',  $
              leq         : '!9l!X',  $ ; `le' is reserved..
              geq         : '!9b!X',  $ ; `ge' is reserved..
              int         : '!MI!X',  $
              partial     : '!MD!X',  $
              nabla       : '!MG!X',  $
              uparrow     : '!97!X',  $
              downarrow   : '!95!X',  $
              leftarrow   : '!94!X',  $
              rightarrow  : '!96!X',  $
              infty       : '!7y!X', $
    ; harder stuff like tilde on top of symbol, etc
    ; Usage: PLOT, [...], TITLE=sym.tilde+'!8Q!X' to get \tilde{Q}, etc
              tilde       : '!S!U!N!A!9A!R!X', $
              bar         : '!S!U!N!A!3-!R!X' $
  }
  endelse

  return, sym

end
