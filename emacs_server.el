;;
;; packages to get
;;

; sudo port install +full

; OR

; sudo port install 
;      texlive-latex
;      texlive-common 
;      texlive-latex-recommended 
;      texlive-latex-extra 
;      ispell
;      ecb 
;      auctex
;      color-theme-mode.el 
;      php-mode.el
;      texlive-fonts-recommended
;      texlive-humanities

;------------------------------------------------------------------------    
; Active
;------------------------------------------------------------------------

;; ;; (setenv "PATH" (getenv "PATH")) ; fix for latex command not working in <MAC> version of Emacs.

(add-to-list 'load-path "~/dropbox/color-theme-6.6.0")
;; ;; (add-to-list 'load-path (expand-file-name "/opt/local/share/emacs/site-lisp/"))
;; ;; (add-to-list 'load-path "~/Emacs/mygoodies/")
;; ;; (add-to-list 'load-path "~/Emacs/matlab")
;; ;; (setenv "PATH" (concat "/opt/local/bin:/usr/local/bin:" (getenv "PATH")))
;; ;; (setq exec-path (append '("/opt/local/bin" "/usr/local/bin") exec-path))


;; ;; (add-to-list 'load-path (expand-file-name "~/Emacs/jdee/lisp"))
;; ;; (add-to-list 'load-path (expand-file-name "~/Emacs/cedet/common"))
;; ;; (load-file (expand-file-name "~/Emacs/cedet/common/cedet.el"))
;; ;; (add-to-list 'load-path (expand-file-name "~/Emacs/elib"))


;------------------------------------------------------------------------
; Inactive
;------------------------------------------------------------------------

;;(add-to-list 'load-path "~/Emacs/tabbar-aquamacs")
;;(add-to-list 'load-path "~/Emacs/ess-5.14/lisp")
;;(add-to-list 'load-path (expand-file-name "~/Emacs/elib"))
;;(add-to-list 'load-path (expand-file-name "~/Emacs/cdet"))


;------------------------------------------------------------------------

;;(server-start)            ; open all files in one instance

(setq ns-pop-up-frames nil) ; Keep files from opening in separate window

;; Set larger font for <MAC>
;;(set-default-font "-apple-Monaco-medium-normal-normal-*-14-*-*-*-m-0-iso10646-1")

;;
;; Drag and Drop on the <MAC>
;;

;;(define-key global-map [ns-drag-file] 'my-ns-open-files)
;;(defun my-ns-open-files ()
;;  "Open files in the list `ns-input-file'."
;;  (interactive)
;;  (mapc 'find-file ns-input-file)
;;  (setq ns-input-file nil))

;;
;; color theme setting
;;

(require 'color-theme)
(eval-after-load "color-theme"
   '(progn
 	(color-theme-initialize)
 	(color-theme-calm-forest)))

;;
;; put scroll bar on the right hand side
;;

;;(set-scroll-bar-mode 'right)

;;(set-scroll-bar-mode 'right)            ; Scroll-bar position: on the
;;(cond ((commandp 'toggle-scroll-bar)    ; right
;;    (toggle-scroll-bar 1)))


(tool-bar-mode -1)               ; Remove top toolbar at the top

(setq inhibit-startup-message t) ; Don't show the startup msg

(fset 'yes-or-no-p 'y-or-n-p)    ; Take y for yes, n for no

;; (require 'tabbar)
;; (tabbar-mode)

;; (require 'linum)                 ; Line Numbers
;; (global-linum-mode 1)

;;
;; Keybindings
;;

(global-set-key "\C-l" 'goto-line)
(global-set-key "\M-[" 'tabbar-backward)
(global-set-key "\M-]" 'tabbar-forward)
(global-set-key "\C-r" 'comment-region)
(global-set-key "\C-t" 'uncomment-region)
(global-set-key "\C-i" 'indent-region)

(global-set-key [M-left] 'windmove-left)          ; move to left windnow
(global-set-key [M-right] 'windmove-right)        ; move to right window
(global-set-key [M-up] 'windmove-up)              ; move to upper window
(global-set-key [M-down] 'windmove-down)          ; move to downer window

;;
;; Turn CUA mode on
;;

(cua-mode t)
    (setq cua-auto-tabify-rectangles nil) ;; Don't tabify after rectangle commands
    (transient-mark-mode 1)               ;; No region when it is not highlighted
    (setq cua-keep-region-after-copy t) 

;;
;; <MAC> key setting
;;
;mac-function-modifier
;mac-control-modifier
;mac-command-modifier
;mac-option-modifier
;mac-right-command-modifier
;mac-right-control-modifier
;mac-right-option-modifier

;;(when (eq system-type 'darwin) ;; mac specific settings
;;  (setq mac-option-modifier 'alt)
;;  (setq mac-command-modifier 'meta)
;;  (global-set-key [kp-delete] 'delete-char) ;; sets fn-delete to be right-delete
;;)

;;(require 'redo)
;;(require 'mac-key-mode)
;;(mac-key-mode 1)

;;
;; copy-paste to terminal
;;

(setq x-select-enable-clipboard t)

;;
;; set tab size
;;

(setq default-tab-width 4)
(define-key global-map (kbd "RET") 'newline-and-indent)

;;
;; split windows vcolcally
;;

(setq display-buffer-prefer-horizontal-split t)
;;(setq split-width-threshold most-positive-fixnum)
;;(setq split-height-threshold 0)

;;
;; save a list of open files in ~/.emacs.desktop
;; save the desktop file automatically if it already exists
;;

(setq desktop-save 'if-exists)
(desktop-save-mode 1)

;; save a bunch of variables to the desktop file
;; for lists specify the len of the maximal saved data also
;; Use M-x desktop-save to keep the opened files persistent between sessions

(setq desktop-globals-to-save
      (append '((extended-command-history . 30)
                (file-name-history        . 100)
                (grep-history             . 30)
                (compile-history          . 30)
                (minibuffer-history       . 50)
                (query-replace-history    . 60)
                (read-expression-history  . 60)
                (regexp-history           . 60)
                (regexp-search-ring       . 20)
                (search-ring              . 20)
                (shell-command-history    . 50)
                tags-file-name
                register-alist)))

;;
;; Enable TRAMP for remote editing over SSH
;;

;; ;; (require 'tramp)
;; ;; (setq tramp-default-method "scp")


;--------------------------------------------------------------------------
; Language Modes
;--------------------------------------------------------------------------


;; MATLAB mode
;; =======================

;; ;; (load-library "matlab-load")
;;(matlab-cedet-setup) ;; Enable CEDET feature support for MATLAB code. (Optional)



;; R mode - "Emacs Speaks Statistics"
;; =======================

;;(require 'ess-site)



;; PhP Mode
;; =======================

;; ;; (load "php-mode")


;; CEDET - Coding Utils
;; =======================

;; ;; (global-ede-mode 1)                      ; Enable the Project management system
;; ;; (semantic-load-enable-code-helpers)      ; Enable prototype help and smart completion 
;; ;; (global-srecode-minor-mode 1)            ; Enable template insertion menu



;; JDEE - Java Env
;; =======================  

;; ;; (require 'jde)



;; ECB - Emacs Code Browser
;; =======================

;;(require 'ecb)
;;(setq ecb-tip-of-the-day nil)
;;(setq ecb-auto-activate t)
;;(require 'ecb-autoloads)



;; AUCTeX mode for (La)TeX
;; =======================

;; ;; (require 'tex-site)


;; AUCTeX specific
;; ===============

;; ;; (load "auctex.el" nil t t)
;; ;; (load "preview-latex.el" nil t t)
;; ;; (setq TeX-parse-self t) ; Enable parse file on load

;;(TeX-auto-save nill)    ; parse info saving
;;(TeX-auto-untabify t)   ; Change tabs to spaces
;;(font-latex-fontify-sectioning 'color)  ; Title highligh. w/ color

;;(defun maximize-frame () 
;;  (interactive)
;;  (set-frame-position (selected-frame) 0 0)
;;  (set-frame-size (selected-frame) 1000 1000))

;;(global-set-key (kbd "<s-S-return>") 'maximize-frame)
