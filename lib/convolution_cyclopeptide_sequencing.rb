#!/usr/bin/ruby

require 'set'

# branch
def expand(peptides, integer_mass)
  new_peptides = Set.new
  peptides.each do |peptide|
    integer_mass.each do |mass|
      # initially peptide is empty set, instead initialise new peptide for each mass
      if peptide.empty? 
        new_peptides << [mass]
      else
        new_peptide = peptide + [mass]
        new_peptides << new_peptide 
      end
    end
  end
  new_peptides
end

# determine mass of peptide by summing over masses
def mass(peptide)
  mass = 0
  for i in 0..peptide.length-1
    mass += peptide[i].to_i
  end
  mass
end

def cyclospectrum(peptide)
  peptide_masses = []
  # always include 0
  whole_peptide = 0
  peptide.each_with_index do |ch, index|
    # sum across each smallest unit for weight of full peptide
    whole_peptide += ch
    # mass of subpeptide at index
    sub = 0
    for i in 0..peptide.length-2
      idx = (index + i) % peptide.length
      sub += peptide[idx]
      peptide_masses << sub
    end
  end
  peptide_masses << whole_peptide
end

def cyclopeptide_scoring(peptide, experimental_spectrum)
  theoretical_spectrum = cyclospectrum(peptide).sort
  #puts theoretical_spectrum.join(' ')
  #puts experimental_spectrum.join(' ')
  t_c = 0
  e_c = 0
  # 0 also counts
  score = 1
  begin
    if experimental_spectrum[e_c] == theoretical_spectrum[t_c]
      score += 1
      e_c += 1
      t_c += 1
    elsif experimental_spectrum[e_c] < theoretical_spectrum[t_c]  
      e_c += 1
    else
      t_c += 1  
    end
  end while !(experimental_spectrum[e_c].nil? || theoretical_spectrum[t_c].nil?)
  score
end

# if peptide is of the form 186-128-113, calculate all linear subpeptides
# e.g. 186+128 = 314 and 128+113 = 241, and 186+128+113 = 427
def cyclospectrum_linear(peptide)
  cyclospectrum = []
  for i in 0..peptide.length-1
    sub = 0
    for j in i..peptide.length-1
      sub += peptide[j]
      cyclospectrum << sub
    end
  end
  cyclospectrum
end

def cyclopeptide_scoring_linear(peptide, experimental_spectrum)
  theoretical_spectrum = cyclospectrum_linear(peptide).sort
  #puts theoretical_spectrum.join(' ')
  #puts experimental_spectrum.join(' ')
  t_c = 0
  e_c = 0
  # 0 also counts
  score = 1
  begin
    #puts "#{score}: #{experimental_spectrum[e_c]} | #{theoretical_spectrum[t_c]}"
    if experimental_spectrum[e_c] == theoretical_spectrum[t_c]
      score += 1
      e_c += 1
      t_c += 1
    elsif experimental_spectrum[e_c] < theoretical_spectrum[t_c]  
      e_c += 1
    else
      t_c += 1  
    end
  end while !(experimental_spectrum[e_c].nil? || theoretical_spectrum[t_c].nil?)
  score
end

# determine at least N highest scoring linear peptides
# ties are included
def trim(leaderboard, spectrum, n)
  lb_map = {}
  leaderboard.each do |peptide|
    linear_score = cyclopeptide_scoring_linear(peptide, spectrum)
    lb_map[peptide] = linear_score
  end
  sorted = lb_map.sort_by {|k,v| v}.reverse
  #puts sorted.join(" ")
  for j in n..leaderboard.size-1
    if sorted[j][1] < sorted[n-1][1]
      #puts "deleting #{sorted[j][0]}"
      leaderboard.delete(sorted[j][0])
    end  
  end
  leaderboard
end

def spectral_convolution(spectrum)
  convolution = []
  for i in 0..spectrum.size-2
    for j in i+1..spectrum.size-1
      diff = spectrum[i] - spectrum[j]
      convolution << (diff).abs unless diff == 0
    end
  end
  #puts convolution.join(' ')
  convolution
end

# select the M most frequent elements between 57 and 200 in the convolution to form an 
# extended alphabet of amino acid masses
def restrict(convolution, m)
  alphabet = []
  freq = convolution.each_with_object({}) { |itm, hsh| hsh[itm] = hsh[itm].to_i + 1 if convolution.include?(itm) }
  eligible = freq.select {|k, v| 57 <= k && k  <= 200}
  sorted = eligible.sort_by {|k,v| v}.reverse
  sorted.each do |am, f|
    if alphabet.size < m
      alphabet << am 
    # let's ensure we include any ties
    elsif alphabet.size == m
      alphabet << am if alphabet[m-1] == am
    end
    break if alphabet[m-1] > am unless alphabet[m-1].nil?
  end
  alphabet
end

# Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.
# Output: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements
# (and ties) of the convolution of Spectrum that fall between 57 and 200, and where the size
# of Leaderboard is restricted to the top N (and ties).
def convolution_cyclopeptide_sequencing(spectrum, n, m, integer_mass)
  convolution = spectral_convolution(spectrum)
  extended_alphabet = restrict(convolution, m)
  #puts extended_alphabet.join(' ')
  leaderboard = Set.new
  leader_peptide = []
  leaderboard << leader_peptide
  leader_score = 0
  mass_parent = spectrum.last.to_i
  puts "mass_parent #{mass_parent}"
  # due to consistency constraints, we can only expand by integer masses that appear
  # in the spectrum and in the extended alphabet
  expandables = extended_alphabet & integer_mass
  #puts expandables.join(' ')
  while !leaderboard.empty?
    leaderboard = expand(leaderboard, expandables)
    leaderboard.each do |peptide|
       mass_peptide = mass(peptide)
      if mass_peptide == mass_parent
        # score the peptide cyclospectrum with experimental spectrum
        pep_score = cyclopeptide_scoring(peptide, spectrum)
        if pep_score > leader_score
          leader_peptide = peptide
          leader_score = pep_score
          puts "new leader with score #{pep_score}: #{leader_peptide.join(' ')}"
        end  
      elsif mass_peptide > mass_parent
        leaderboard.delete(peptide)
      end
    end
    leaderboard = trim(leaderboard.to_set, spectrum, n)
  end
  leader_peptide
end

m = 16
n = 372
integer_mass = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
spectrum = "57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493".split(" ").to_a.collect{|s| s.to_i}.sort

leaderpeptide = convolution_cyclopeptide_sequencing(spectrum, n, m, integer_mass)
puts leaderpeptide.join('-')