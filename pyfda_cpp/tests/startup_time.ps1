# Start time of pyfda_cpp.exe and pyfdax.exe (PyInstaller, PR #4) on Windows: time from the
# start of the process until the main window exists and until it is ready for input
# (its event loop answers a message). The first run after the download is "cold", the following ones "warm".
param([string]$Cpp, [string]$Py, [int]$Runs = 10)

Add-Type @"
using System;
using System.Runtime.InteropServices;
public static class W32 {
    [DllImport("user32.dll", CharSet = CharSet.Unicode)] public static extern IntPtr FindWindow(string cls, string title);
    [DllImport("user32.dll")] public static extern bool IsWindowVisible(IntPtr h);
    [DllImport("user32.dll")] public static extern IntPtr SendMessageTimeout(IntPtr h, uint msg, IntPtr w, IntPtr l, uint flags, uint ms, out IntPtr res);
}
"@

# the main window is found by its title (pyfdax.exe also has a console window)
function Measure-Start([string]$exe, [string]$name, [string]$title, [string[]]$argv) {
    $sw = [Diagnostics.Stopwatch]::StartNew()
    if ($argv) { $p = Start-Process -FilePath $exe -ArgumentList $argv -PassThru } else { $p = Start-Process -FilePath $exe -PassThru }
    $h = [IntPtr]::Zero
    while ($sw.Elapsed.TotalSeconds -lt 120) {
        $h = [W32]::FindWindow($null, $title)
        if ($h -ne [IntPtr]::Zero -and [W32]::IsWindowVisible($h)) { break }
        Start-Sleep -Milliseconds 5
    }
    $t_win = $sw.Elapsed.TotalSeconds
    $t_ready = [double]::NaN
    $found = $h -ne [IntPtr]::Zero
    if ($found) {
        # WM_NULL returns once the event loop of the window processes messages
        $r = [IntPtr]::Zero
        [void][W32]::SendMessageTimeout($h, 0, [IntPtr]::Zero, [IntPtr]::Zero, 0, 60000, [ref]$r)
        $t_ready = $sw.Elapsed.TotalSeconds
    }
    Get-Process -Name $name -ErrorAction SilentlyContinue | Stop-Process -Force
    Start-Sleep -Milliseconds 1000
    return [pscustomobject]@{ window = $t_win; ready = $t_ready; found = $found }
}

$cfg = Join-Path $env:RUNNER_TEMP 'pyfda_cpp_cfg'
$rows = @()
foreach ($app in @(@{ n = 'pyfda_cpp'; exe = $Cpp; t = 'pyfda C++ - Python Filter Design Analysis Tool'; a = @('--config-dir', $cfg) },
                   @{ n = 'pyfdax'; exe = $Py; t = 'pyfda - Python Filter Design and Analysis'; a = $null })) {
    for ($i = 0; $i -le $Runs; $i++) {
        $r = Measure-Start $app.exe $app.n $app.t $app.a
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
