.text
.globl readSignalFromFile
.globl writeOutputToFile
.globl applyWienerFilter
.globl computeMMSE

# FUNCTION: readSignalFromFile
# DESCRIPTION: Opens the "input.txt" file, reads its entire content into a 
# text buffer, and calls the parsing subroutine to extract the data into 
# the 'desired_signal' and 'input_signal' arrays.
readSignalFromFile:
    addi    $sp, $sp, -4
    sw      $ra, 0($sp)

    li      $v0, 13
    la      $a0, input_file
    li      $a1, 0
    li      $a2, 0
    syscall
    move    $s0, $v0

    li      $v0, 14
    move    $a0, $s0
    la      $a1, buffer_read
    li      $a2, 2048
    syscall
    
    li      $v0, 16
    move    $a0, $s0
    syscall

    la      $t0, buffer_read
    la      $a1, desired_signal
    jal     parseStringToArray
    beqz    $v0, end_read_file

    la      $a1, input_signal
    jal     parseStringToArray

end_read_file:
    lw      $ra, 0($sp)
    addi    $sp, $sp, 4
    jr      $ra

# FUNCTION: parseStringToArray
# DESCRIPTION: Iterates through the text buffer to identify numbers. It handles 
# negative signs, integer parts, and decimal fractions, converting the ASCII 
# characters into IEEE 754 floating-point values and storing them in an array.
parseStringToArray:
    li      $t1, 0
    li      $t2, 10
    lwc1    $f20, float_ten
    
parse_loop:
    bge     $t1, $t2, parse_success
    li      $t3, 0      
    lwc1    $f0, float_zero
    lwc1    $f2, float_one
    
skip_space:
    lb      $t4, 0($t0)
    addi    $t0, $t0, 1
    beq     $t4, 32, skip_space
    beq     $t4, 10, skip_space
    beq     $t4, 13, skip_space
    beqz    $t4, parse_fail
    
    bne     $t4, 45, check_digit
    li      $t3, 1
    lb      $t4, 0($t0)
    addi    $t0, $t0, 1

check_digit:
    beq     $t4, 46, parse_fraction
    beq     $t4, 32, store_val
    beq     $t4, 10, store_val
    beq     $t4, 13, store_val
    beqz    $t4, store_val
    
    sub     $t4, $t4, 48
    mtc1    $t4, $f4
    cvt.s.w $f4, $f4
    
    mul.s   $f0, $f0, $f20
    add.s   $f0, $f0, $f4
    
    lb      $t4, 0($t0)
    addi    $t0, $t0, 1
    j       check_digit

parse_fraction:
    lb      $t4, 0($t0)
    addi    $t0, $t0, 1
    beq     $t4, 32, store_val
    beq     $t4, 10, store_val
    beq     $t4, 13, store_val
    beqz    $t4, store_val
    
    sub     $t4, $t4, 48
    mtc1    $t4, $f4
    cvt.s.w $f4, $f4
    
    mul.s   $f2, $f2, $f20
    div.s   $f4, $f4, $f2
    add.s   $f0, $f0, $f4
    j       parse_fraction

store_val:
    beqz    $t3, save_to_array
    neg.s   $f0, $f0

save_to_array:
    sll     $t5, $t1, 2
    add     $t5, $a1, $t5
    swc1    $f0, 0($t5)
    addi    $t1, $t1, 1
    j       parse_loop

parse_fail:
    li      $v0, 0
    jr      $ra

parse_success:
    li      $v0, 1
    jr      $ra

# FUNCTION: applyWienerFilter
# DESCRIPTION: Applies the calculated Wiener filter coefficients to the input 
# signal using convolution. It computes y(n) = sum(h[k] * x[n-k]) for M=3.
applyWienerFilter:
    li      $t0, 0
    li      $t3, 3
    li      $t4, 10
    
loop_n_filter:
    bge     $t0, $t4, end_filter
    lwc1    $f0, float_zero
    li      $t1, 0

loop_k_filter:
    bge     $t1, $t3, end_k_filter
    sub     $t2, $t0, $t1
    bltz    $t2, next_k_filter
    
    sll     $t5, $t1, 2
    add     $t5, $a1, $t5
    lwc1    $f2, 0($t5)
    
    sll     $t6, $t2, 2
    add     $t6, $a0, $t6
    lwc1    $f4, 0($t6)
    
    mul.s   $f6, $f2, $f4
    add.s   $f0, $f0, $f6

next_k_filter:
    addi    $t1, $t1, 1
    j       loop_k_filter

end_k_filter:
    sll     $t7, $t0, 2
    add     $t7, $a2, $t7
    swc1    $f0, 0($t7)
    addi    $t0, $t0, 1
    j       loop_n_filter

end_filter:
    jr      $ra

# FUNCTION: computeMMSE
# DESCRIPTION: Calculates the Minimum Mean Square Error (MMSE). It first 
# computes the variance of the desired signal, then subtracts the product 
# of the optimal coefficients and the cross-correlation values.
computeMMSE:
    la      $t0, desired_signal
    li      $t1, 0
    li      $t4, 10
    lwc1    $f0, float_zero
    
loop_var_d:
    bge     $t1, $t4, calc_var_d
    sll     $t5, $t1, 2
    add     $t5, $t0, $t5
    lwc1    $f2, 0($t5)
    mul.s   $f4, $f2, $f2
    add.s   $f0, $f0, $f4
    addi    $t1, $t1, 1
    j       loop_var_d

calc_var_d:
    lwc1    $f10, float_ten
    div.s   $f0, $f0, $f10
    
    la      $t2, optimize_coefficient
    la      $t3, rdx_array
    li      $t1, 0
    li      $t4, 3

loop_mmse:
    bge     $t1, $t4, end_mmse
    sll     $t5, $t1, 2
    add     $t6, $t2, $t5
    lwc1    $f6, 0($t6)
    add     $t7, $t3, $t5
    lwc1    $f8, 0($t7)
    mul.s   $f12, $f6, $f8
    sub.s   $f0, $f0, $f12
    addi    $t1, $t1, 1
    j       loop_mmse

end_mmse:
    abs.s   $f0, $f0
    la      $t8, mmse
    swc1    $f0, 0($t8)
    jr      $ra

# FUNCTION: writeOutputToFile
# DESCRIPTION: Prints the filtered output array and MMSE value to the console. 
# It also builds an ASCII buffer containing the same formatted information 
# and writes it into the "output.txt" file.
writeOutputToFile:
    addi    $sp, $sp, -4
    sw      $ra, 0($sp)

    li      $v0, 4
    la      $a0, msg_filtered
    syscall
    
    la      $s1, output_signal
    li      $s2, 0
    li      $s3, 10
    la      $s4, buffer_write
    
    la      $t0, msg_filtered
copy_msg1:
    lb      $t1, 0($t0)
    beqz    $t1, loop_print_out
    sb      $t1, 0($s4)
    addi    $t0, $t0, 1
    addi    $s4, $s4, 1
    j       copy_msg1

loop_print_out:
    bge     $s2, $s3, print_mmse_label
    sll     $t5, $s2, 2
    add     $t5, $s1, $t5
    lwc1    $f12, 0($t5)
    
    li      $v0, 2
    syscall
    li      $v0, 11
    li      $a0, 32
    syscall

    jal     float_to_string
    li      $t1, 32
    sb      $t1, 0($s4)
    addi    $s4, $s4, 1
    
    addi    $s2, $s2, 1
    j       loop_print_out

print_mmse_label:
    li      $v0, 4
    la      $a0, msg_mmse
    syscall
    
    la      $t8, mmse
    lwc1    $f12, 0($t8)
    li      $v0, 2
    syscall
    
    li      $v0, 11
    li      $a0, 10
    syscall

    la      $t0, msg_mmse
copy_msg2:
    lb      $t1, 0($t0)
    beqz    $t1, print_mmse_val
    sb      $t1, 0($s4)
    addi    $t0, $t0, 1
    addi    $s4, $s4, 1
    j       copy_msg2

print_mmse_val:
    la      $t8, mmse
    lwc1    $f12, 0($t8)
    jal     float_to_string

    li      $v0, 13
    la      $a0, output_file
    li      $a1, 1
    li      $a2, 0
    syscall
    move    $s0, $v0
    
    la      $t0, buffer_write
    sub     $a2, $s4, $t0
    li      $v0, 15
    move    $a0, $s0
    la      $a1, buffer_write
    syscall
    
    li      $v0, 16
    move    $a0, $s0
    syscall

    lw      $ra, 0($sp)
    addi    $sp, $sp, 4
    jr      $ra

# FUNCTION: float_to_string (Helper)
# DESCRIPTION: Converts a floating-point number into an ASCII string with one 
# decimal place and appends it to the global text buffer.
float_to_string:
    lwc1    $f0, float_zero
    c.lt.s  $f12, $f0
    bc1f    f2s_positive
    li      $t1, 45
    sb      $t1, 0($s4)
    addi    $s4, $s4, 1
    abs.s   $f12, $f12
    
f2s_positive:
    lwc1    $f0, float_ten
    mul.s   $f12, $f12, $f0
    round.w.s $f12, $f12
    mfc1    $t0, $f12
    
    li      $t1, 10
    div     $t0, $t1
    mflo    $t2
    mfhi    $t3
    
    li      $t4, 0
f2s_int_loop:
    div     $t2, $t1
    mflo    $t2
    mfhi    $t5
    addi    $t5, $t5, 48
    addi    $sp, $sp, -4
    sw      $t5, 0($sp)
    addi    $t4, $t4, 1
    bnez    $t2, f2s_int_loop
    
f2s_pop_loop:
    lw      $t5, 0($sp)
    addi    $sp, $sp, 4
    sb      $t5, 0($s4)
    addi    $s4, $s4, 1
    addi    $t4, $t4, -1
    bnez    $t4, f2s_pop_loop
    
    li      $t5, 46
    sb      $t5, 0($s4)
    addi    $s4, $s4, 1
    
    addi    $t3, $t3, 48
    sb      $t3, 0($s4)
    addi    $s4, $s4, 1
    jr      $ra

# DATA SECTION
# DESCRIPTION: Declares strings, file names, floating-point constants, and 
# memory spaces required for reading/writing files and signal processing.
.data
    buffer_read:          .space 2048
    buffer_write:         .space 1024
    float_zero:           .float 0.0
    float_one:            .float 1.0
    float_ten:            .float 10.0
    input_file:           .asciiz "input.txt"
    output_file:          .asciiz "output.txt"
    msg_filtered:         .asciiz "Filtered output: "
    msg_mmse:             .asciiz "\nMMSE: "
    
    # Global Arrays (Integration placeholders)
    desired_signal:       .space 40
    input_signal:         .space 40
    output_signal:        .space 40
    optimize_coefficient: .space 12
    rdx_array:            .space 12
    mmse:                 .float 0.0