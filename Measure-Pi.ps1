[CmdletBinding()]
Param (
    [Parameter(Mandatory=$true)] [string] $Terms,
    [Parameter(Mandatory=$true)] [string] $Method,
	[Parameter(Mandatory=$true)] [string] $Path,
	[Parameter(Mandatory=$true)] [string] $NumberOfProcess
)

Measure-Command {
	Write-Host '3.141592653589793238462643(Aim)'
	Write-Host '3.1415926535897932384626433832795028841971693993751(50 digits)'
	mpiexec -n $NumberOfProcess $Path -t $Terms -m $Method | Write-Host
	Write-Host '          11111111112222222222333333333344444444445'
	Write-Host '1 2345678901234567890123456789012345678901234567890'
}
