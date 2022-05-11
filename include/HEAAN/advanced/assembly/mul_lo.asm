; ----------------------------------------------------------------------------------------------------------- 
; Text or code segment
; ----------------------------------------------------------------------------------------------------------- 

_TEXT SEGMENT
  
; The PUBLIC modifier will make your function visible and callable outside
PUBLIC mul_lo_asm
 
mul_lo_asm PROC

; rcx : a, rdx : b
mov rax, rcx
mul rdx

 ret 
mul_lo_asm ENDP
 
_TEXT ENDS
 
END