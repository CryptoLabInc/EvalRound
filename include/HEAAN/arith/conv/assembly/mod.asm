; ----------------------------------------------------------------------------------------------------------- 
; Text or code segment
; ----------------------------------------------------------------------------------------------------------- 

_TEXT SEGMENT
  
; The PUBLIC modifier will make your function visible and callable outside
PUBLIC mod_asm
 
mod_asm PROC

; rcx : al, rdx : ah, r8 : q

; rax <- al
; rcx <- ah
; rdx <- 0
; rdx <- al % q
; r9 <- rdx
mov rax, rcx
mov rcx, rdx
xor rdx, rdx
div r8
mov r9, rdx

; rax <- ah
; rdx <- 0
; rdx <- ah % q
mov rax, rcx
xor rdx, rdx
div r8

; rax <- al % q
; rdx <- ah % q
; rdx <- res
; rax <- rdx
mov rax, r9
div r8
mov rax, rdx
 
 ret 
mod_asm ENDP
 
_TEXT ENDS
 
END