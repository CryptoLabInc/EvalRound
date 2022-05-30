; ----------------------------------------------------------------------------------------------------------- 
; Text or code segment
; ----------------------------------------------------------------------------------------------------------- 

_TEXT SEGMENT
  
; The PUBLIC modifier will make your function visible and callable outside
PUBLIC mul_hi_asm
 
mul_hi_asm PROC

; rcx : a, rdx : b
mov rax, rcx
mul rdx
mov rax, rdx

 ret 
mul_hi_asm ENDP
 
_TEXT ENDS
 
END