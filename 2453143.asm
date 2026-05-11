# =====================================================
#  Le Tran Xuan Tao - 2453143
#  Description:
#     - computeAutocorrelation
#     - computeCrosscorrelation
#     - createToeplitzMatrix
#  Float precision, N = 3
# =====================================================

.data                               # bat dau vùng du lieu
N:           .word 3                # khai báo N = 3

# --- Input signals ---
signal:      .float 1.0, 2.0, 3.0  # mang tín hieu input x[n]
desired:     .float 1.5, 2.8, 3.2  # mang tín hieu mong muon d[n]

# --- Output buffers ---
autocorr:    .space 12              # cap phát 12 byte cho autocorrelation (3 float)
crosscorr:   .space 12              # cap phát 12 byte cho crosscorrelation (3 float)
R:           .space 36              # cap phát 36 byte cho ma tran Toeplitz 3x3

# --- UI strings ---
msg_auto:    .asciiz "\n--- Autocorrelation ---\n" # chu?i tiêu ?? autocorrelation
msg_cross:   .asciiz "\n--- Crosscorrelation ---\n" # chu?i tiêu ?? crosscorrelation
msg_toep:    .asciiz "\n--- Toeplitz Matrix ---\n" # chu?i tiêu ?? ma tr?n Toeplitz
newline:     .asciiz "\n"          # ký t? xu?ng dòng

.text                               # b?t ??u vùng code
j main                              # nh?y ??n hàm main
.globl main                         # khai báo main là global

# =====================================================
# computeAutocorrelation
# =====================================================
computeAutocorrelation:             # b?t ??u hàm computeAutocorrelation
#T?o stack frame, l?u c?c thanh ghi t?m
    addi $sp, $sp, -12              # c?p phát 12 byte stack
    sw   $ra, 8($sp)                # l?u ??a ch? tr? v?
    sw   $s0, 4($sp)                # l?u thanh ghi s0
    sw   $s1, 0($sp)                # l?u thanh ghi s1

    move $t7, $a2                   # copy N t? a2 sang t7
    li   $s0, 0                     # kh?i t?o k = 0

k_loop:                             # vòng l?p ngoài theo k
    bge  $s0, $t7, autocorr_done    # n?u k >= N thì k?t thúc
    mtc1 $zero, $f0                 # chuy?n s? nguyên 0 sang thanh ghi float
    cvt.s.w $f0, $f0                # ép ki?u int sang float ?? sum = 0.0
    move $s1, $s0                   # gán i = k

i_loop:                             # vòng l?p trong theo i
    bge  $s1, $t7, i_done           # n?u i >= N thì k?t thúc vòng l?p trong
    sll  $t0, $s1, 2                # tính offset i*4
    add  $t0, $t0, $a0              # c?ng base address signal
    lwc1 $f2, 0($t0)                # load x(i)

    sub  $t1, $s1, $s0              # tính i-k
    sll  $t1, $t1, 2                # tính offset (i-k)*4
    add  $t1, $t1, $a0              # c?ng base address signal
    lwc1 $f4, 0($t1)                # load x(i-k)

    mul.s $f6, $f2, $f4             # tính x(i)*x(i-k)
    add.s $f0, $f0, $f6             # c?ng vào sum

    addi $s1, $s1, 1                # i++
    j i_loop                        # quay l?i ??u vòng l?p trong

i_done:                             # k?t thúc vòng l?p trong
    mtc1 $t7, $f8                   # chuy?n N sang thanh ghi float
    cvt.s.w $f8, $f8                # ép ki?u N thành float
    div.s $f0, $f0, $f8             # chia sum cho N

    sll  $t3, $s0, 2                # tính offset k*4
    add  $t3, $t3, $a1              # c?ng base address autocorr
    swc1 $f0, 0($t3)                # l?u autocorr[k]

    addi $s0, $s0, 1                # k++
    j k_loop                        # quay l?i vòng l?p ngoài

autocorr_done:                      # k?t thúc hàm autocorrelation
#ph?c h?i tr?ng thái ban ??u
    lw   $ra, 8($sp)                # khôi ph?c ??a ch? tr? v?
    lw   $s0, 4($sp)                # khôi ph?c s0
    lw   $s1, 0($sp)                # khôi ph?c s1
    addi $sp, $sp, 12               # gi?i phóng stack
    jr   $ra                        # quay v? hàm g?i

# =====================================================
# computeCrosscorrelation (theo  Wiener)
# rxd(l) = (1/N) * sum_{n=0}^{N-l-1} x(n+l)*d(n)
# =====================================================
computeCrosscorrelation:            # b?t ??u hàm computeCrosscorrelation
    addi $sp, $sp, -20              # c?p phát stack frame
    sw   $ra, 16($sp)               # l?u ??a ch? tr? v?
    sw   $s0, 12($sp)               # l?u s0
    sw   $s1, 8($sp)                # l?u s1
    sw   $s2, 4($sp)                # l?u s2
    sw   $s3, 0($sp)                # l?u s3

    move $s2, $a3                   # copy N t? a3 sang s2
    li   $s0, 0                     # kh?i t?o l = 0

outer_l:                            # vòng l?p ngoài theo l
    bge  $s0, $s2, cc_end           # n?u l >= N thì k?t thúc
    mtc1 $zero, $f0                 # gán 0 vào float register
    cvt.s.w $f0, $f0                # sum = 0.0

    li   $s1, 0                     # kh?i t?o n = 0
inner_n:                            # vòng l?p trong theo n
    sub  $t0, $s2, $s0              # tính N-l
    bge  $s1, $t0, sum_done         # n?u n >= N-l thì k?t thúc

    add  $t1, $s1, $s0              # tính n+l
    sll  $t1, $t1, 2                # offset (n+l)*4
    add  $t1, $t1, $a1              # c?ng base address x[]
    lwc1 $f2, 0($t1)                # load x(n+l)

    sll  $t2, $s1, 2                # offset n*4
    add  $t2, $t2, $a0              # c?ng base address d[]
    lwc1 $f4, 0($t2)                # load d(n)

    mul.s $f6, $f2, $f4             # nhân x(n+l)*d(n)
    add.s $f0, $f0, $f6             # c?ng vào sum

    addi $s1, $s1, 1                # n++
    j inner_n                       # quay l?i vòng l?p trong

sum_done:                           # k?t thúc tính t?ng
    mtc1 $s2, $f8                   # chuy?n N sang float register
    cvt.s.w $f8, $f8                # ép ki?u float
    div.s $f0, $f0, $f8             # chia sum cho N

    sll  $t4, $s0, 2                # tính offset l*4
    add  $t4, $t4, $a2              # c?ng base address crosscorr
    swc1 $f0, 0($t4)                # l?u crosscorr[l]

    addi $s0, $s0, 1                # l++
    j outer_l                       # quay l?i vòng l?p ngoài

cc_end:                             # k?t thúc hàm crosscorrelation
    lw   $ra, 16($sp)               # khôi ph?c ??a ch? tr? v?
    lw   $s0, 12($sp)               # khôi ph?c s0
    lw   $s1, 8($sp)                # khôi ph?c s1
    lw   $s2, 4($sp)                # khôi ph?c s2
    lw   $s3, 0($sp)                # khôi ph?c s3
    addi $sp, $sp, 20               # gi?i phóng stack
    jr   $ra                        # quay v? hàm g?i

# =====================================================
# createToeplitzMatrix
# =====================================================
createToeplitzMatrix:               # b?t ??u hàm createToeplitzMatrix
#T?o stack frame
    addi $sp, $sp, -16              # c?p phát stack frame
    sw   $ra, 12($sp)               # l?u ??a ch? tr? v?
    sw   $s0, 8($sp)                # l?u s0
    sw   $s1, 4($sp)                # l?u s1
    sw   $s2, 0($sp)                # l?u s2

    move $s2, $a2                   # copy N sang s2
    li   $s0, 0                     # kh?i t?o i = 0

outer_loop:                         # vòng l?p ngoài theo i
    bge  $s0, $s2, end_outer        # n?u i >= N thì k?t thúc
    li   $s1, 0                     # kh?i t?o j = 0

inner_loop:                         # vòng l?p trong theo j
    bge  $s1, $s2, end_inner        # n?u j >= N thì k?t thúc vòng l?p trong

    sub  $t0, $s0, $s1              # tính i-j
    bgez $t0, abs_ok                # n?u >=0 thì b? qua
    sub  $t0, $zero, $t0            # l?y tr? tuy?t ??i n?u âm
abs_ok:                             # nhãn x? lý tr? tuy?t ??i

    blt  $t0, $s2, within_bounds    # n?u |i-j| < N thì h?p l?
    addi $t0, $s2, -1               # n?u v??t biên thì gán = N-1

within_bounds:                      # nhãn x? lý h?p l?
    sll  $t1, $t0, 2                # tính offset |i-j|*4
    add  $t1, $t1, $a0              # c?ng base autocorr
    lwc1 $f2, 0($t1)                # load autocorr[|i-j|]

    mul  $t2, $s0, $s2              # tính i*N
    add  $t2, $t2, $s1              # tính i*N+j
    sll  $t2, $t2, 2                # ??i sang byte offset
    add  $t2, $t2, $a1              # c?ng base address R
    swc1 $f2, 0($t2)                # l?u vào R[i][j]

    addi $s1, $s1, 1                # j++
    j inner_loop                    # quay l?i vòng l?p trong

end_inner:                          # k?t thúc vòng l?p trong
    addi $s0, $s0, 1                # i++
    j outer_loop                    # quay l?i vòng l?p ngoài

end_outer:                          # k?t thúc hàm Toeplitz
#Ph?c h?i
    lw   $ra, 12($sp)               # khôi ph?c ??a ch? tr? v?
    lw   $s0, 8($sp)                # khôi ph?c s0
    lw   $s1, 4($sp)                # khôi ph?c s1
    lw   $s2, 0($sp)                # khôi ph?c s2
    addi $sp, $sp, 16               # gi?i phóng stack
    jr $ra                          # quay v? hàm g?i

# =====================================================
# MAIN: ch?y full pipeline
# =====================================================
main:                               # b?t ??u hàm main
    lw   $a2, N                     # load N vào a2
    la   $a0, signal                # load ??a ch? signal
    la   $a1, autocorr              # load ??a ch? autocorr
    jal  computeAutocorrelation     # g?i hàm autocorrelation

    lw   $a3, N                     # load N vào a3
    la   $a0, desired               # load ??a ch? desired
    la   $a1, signal                # load ??a ch? signal
    la   $a2, crosscorr             # load ??a ch? crosscorr
    jal  computeCrosscorrelation    # g?i hàm crosscorrelation

    lw   $a2, N                     # load N
    la   $a0, autocorr              # load ??a ch? autocorr
    la   $a1, R                     # load ??a ch? ma tr?n R
    jal  createToeplitzMatrix       # g?i hàm t?o Toeplitz

    # --- Print autocorr ---
    li   $v0, 4                     # syscall print string
    la   $a0, msg_auto              # load chu?i tiêu ??
    syscall                         # in tiêu ??

    li $t0, 0                       # kh?i t?o bi?n ??m
    lw $t1, N                       # load N
print_auto:                         # vòng l?p in autocorr
    bge $t0, $t1, print_cross       # n?u h?t m?ng thì chuy?n sang crosscorr
    sll $t2, $t0, 2                 # tính offset
    la  $t3, autocorr               # load base autocorr
    add $t3, $t3, $t2               # c?ng offset
    lwc1 $f12, 0($t3)               # load float c?n in
    li   $v0, 2                     # syscall print float
    syscall                         # in float
    li   $v0, 4                     # syscall print string
    la   $a0, newline               # load newline
    syscall                         # in newline
    addi $t0, $t0, 1                # t?ng bi?n ??m
    j print_auto                    # quay l?i vòng l?p

# --- Print crosscorr ---
print_cross:                        # nhãn in crosscorr
    li   $v0, 4                     # syscall print string
    la   $a0, msg_cross             # load tiêu ??
    syscall                         # in tiêu ??

    li $t0, 0                       # reset bi?n ??m
    lw $t1, N                       # load N
print_cross_loop:                   # vòng l?p in crosscorr
    bge $t0, $t1, print_toep        # n?u h?t thì chuy?n sang Toeplitz
    sll $t2, $t0, 2                 # tính offset
    la  $t3, crosscorr              # load base crosscorr
    add $t3, $t3, $t2               # c?ng offset
    lwc1 $f12, 0($t3)               # load float
    li   $v0, 2                     # syscall print float
    syscall                         # in float
    li   $v0, 4                     # syscall print string
    la   $a0, newline               # load newline
    syscall                         # in newline
    addi $t0, $t0, 1                # t?ng bi?n ??m
    j print_cross_loop              # l?p l?i

# --- Print Toeplitz matrix ---
print_toep:                         # nhãn in Toeplitz
    li   $v0, 4                     # syscall print string
    la   $a0, msg_toep              # load tiêu ??
    syscall                         # in tiêu ??

    lw $t6, N                       # load N
    li $t0, 0                       # i = 0
outer_print:                        # vòng l?p theo hàng
    bge $t0, $t6, done              # n?u i >= N thì k?t thúc
    li $t1, 0                       # j = 0
inner_print:                        # vòng l?p theo c?t
    bge $t1, $t6, next_row          # n?u j >= N thì sang hàng m?i
    mul $t2, $t0, $t6               # tính i*N
    add $t2, $t2, $t1               # tính i*N+j
    sll $t2, $t2, 2                 # ??i sang byte offset
    la  $t3, R                      # load base ma tr?n R
    add $t3, $t3, $t2               # c?ng offset
    lwc1 $f12, 0($t3)               # load ph?n t? R[i][j]
    li   $v0, 2                     # syscall print float
    syscall                         # in float

    li   $v0, 11                    # syscall print char
    li   $a0, ' '                   # ký t? kho?ng tr?ng
    syscall                         # in kho?ng tr?ng

    addi $t1, $t1, 1                # j++
    j inner_print                   # quay l?i vòng l?p c?t

next_row:                           # sang hàng ti?p theo
    li   $v0, 4                     # syscall print string
    la   $a0, newline               # load newline
    syscall                         # in newline
    addi $t0, $t0, 1                # i++
    j outer_print                   # quay l?i vòng l?p hàng

done:                               # k?t thúc ch??ng trình
    li $v0, 10                      # syscall exit
    syscall                         # thoát ch??ng trình
