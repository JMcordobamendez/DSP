# Start time of pyfda_cpp.exe and pyfdax.exe (PyInstaller, PR #4) on Windows: time from the
# start of the process until the main window exists and until it is ready for input
# (WaitForInputIdle). The first run after the download is "cold", the following ones "warm".
param([string]$Cpp, [string]$Py, [int]$Runs = 10)

function Measure-Start([string]$exe, [string]$name, [string[]]$argv) {
    $sw = [Diagnostics.Stopwatch]::StartNew()
    if ($argv) { $p = Start-Process -FilePath $exe -ArgumentList $argv -PassThru } else { $p = Start-Process -FilePath $exe -PassThru }
    $win = $null
    while ($sw.Elapsed.TotalSeconds -lt 120) {
        # PyInstaller onefile: the window belongs to a child process with the same name
        $win = Get-Process -Name $name -ErrorAction SilentlyContinue | Where-Object { $_.MainWindowHandle -ne 0 } | Select-Object -First 1
        if ($win) { break }
        Start-Sleep -Milliseconds 10
    }
    $t_win = $sw.Elapsed.TotalSeconds
    $t_idle = [double]::NaN
    if ($win) { [void]$win.WaitForInputIdle(60000); $t_idle = $sw.Elapsed.TotalSeconds }
    Get-Process -Name $name -ErrorAction SilentlyContinue | Stop-Process -Force
    Start-Sleep -Milliseconds 500
    return [pscustomobject]@{ window = $t_win; ready = $t_idle; found = [bool]$win }
}

$cfg = Join-Path $env:RUNNER_TEMP 'pyfda_cpp_cfg'
$rows = @()
foreach ($app in @(@{ n = 'pyfda_cpp'; exe = $Cpp; a = @('--config-dir', $cfg) }, @{ n = 'pyfdax'; exe = $Py; a = $null })) {
    for ($i = 0; $i -le $Runs; $i++) {
        $r = Measure-Start $app.exe $app.n $app.a
        $rows += [pscustomobject]@{ app = $app.n; run = $(if ($i -eq 0) { 'cold' } else { "warm $i" }); window_s = [math]::Round($r.window, 3); ready_s = [math]::Round($r.ready, 3); found = $r.found }
    }
}
$rows | Format-Table -AutoSize | Out-String -Width 200 | Write-Output
foreach ($n in 'pyfda_cpp', 'pyfdax') {
    $w = $rows | Where-Object { $_.app -eq $n -and $_.run -ne 'cold' } | Sort-Object ready_s
    $c = $rows | Where-Object { $_.app -eq $n -and $_.run -eq 'cold' }
    $med = $w[[int][math]::Floor($w.Count / 2)].ready_s
    Write-Output ("{0}: cold {1} s, warm median {2} s (min {3}, max {4})" -f $n, $c.ready_s, $med, $w[0].ready_s, $w[-1].ready_s)
}
$rows | Export-Csv -NoTypeInformation startup_time.csv
