subroutine cgns_test()

#IFDEF USE_CGNS 

print *,'CGNS enabled in this build'
  
#ELSE

print *,'CGNS not enabled in this build'

#ENDIF

end subroutine cgns_test
