subroutine cgns_test()

#IFDEF USE_CGNS 


  
#ELSE

print *,'CGNS not enabled in this build'

#ENDIF

end subroutine cgns_test
